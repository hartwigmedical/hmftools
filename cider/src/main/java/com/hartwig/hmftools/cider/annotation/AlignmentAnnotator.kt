package com.hartwig.hmftools.cider.annotation

import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.IgTcrGene.Companion.fromCommonIgTcrGene
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_MATCH_REF_IDENTITY
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_ALIGN_SCORE_MIN
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_VDJ_FLANK_BASES
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_VJ_IDENTITY_MIN
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
import com.hartwig.hmftools.common.cider.IgTcrRegion
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.CigarElement
import htsjdk.samtools.util.SequenceUtil.reverseComplement
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.math.max
import kotlin.math.min

enum class AlignmentStatus
{
    V_D_J, V_J, V_D, D_J, V_ONLY, D_ONLY, J_ONLY, NO_REARRANGEMENT, NO_VDJ_ALIGNMENT, SKIPPED_ALIGN
}

data class GeneAnnotation(
    val gene: IgTcrGene,
    val alignment: Alignment,
    val identity: Double?,
    val supplementaryGenes: List<IgTcrGene>
)

data class AlignmentAnnotation(
    val vdjSequence: VDJSequence,
    val status: AlignmentStatus,
    val fullAlignment: Alignment? = null,
    val vGene: GeneAnnotation? = null,
    val dGene: GeneAnnotation? = null,
    val jGene: GeneAnnotation? = null
)

class AlignmentAnnotator
{
    private val sLogger = LogManager.getLogger(AlignmentAnnotator::class.java)

    private val mRefGenomeVersion: RefGenomeVersion
    private val mRefGenomeDictPath: String
    private val mRefGenomeBwaIndexImagePath: String
    private val mVdjGenes: Map<String, IgTcrGene>
    private val mImgtSequences: ImgtSequenceFile

    constructor(refGenomeVersion: RefGenomeVersion, refGenomeDictPath: String, refGenomeBwaIndexImagePath: String)
    {
        mRefGenomeVersion = refGenomeVersion
        mRefGenomeDictPath = refGenomeDictPath
        mRefGenomeBwaIndexImagePath = refGenomeBwaIndexImagePath

        val igTcrGenes = IgTcrGeneFile.read(refGenomeVersion)
            .map(::fromCommonIgTcrGene)
            .filter { it.region.isVDJ }
        val vdjGenes: HashMap<String, IgTcrGene> = HashMap()
        for (geneData in igTcrGenes)
        {
            val key = geneData.geneAllele
            if (key in vdjGenes)
            {
                throw RuntimeException("Duplicate gene allele: $key")
            }
            vdjGenes[key] = geneData
        }
        mVdjGenes = vdjGenes

        mImgtSequences = ImgtSequenceFile(refGenomeVersion)
    }

    fun runAnnotate(sampleId: String, vdjList: List<VDJSequence>, outputDir: String, numThreads: Int)
            : Collection<AlignmentAnnotation>
    {
        sLogger.info("Running alignment annotation")

        val alignmentMetadata = vdjList.map {
            val querySeqRange = alignmentQuerySeqRange(it)
            AlignmentMetadata(it, querySeqRange, it.layout.consensusSequenceString().substring(querySeqRange))
        }
        val querySequences = alignmentMetadata.map { it.querySeq }

        // Align to the normal reference genome to get most of the alignments.
        // For GRCh37, also align to a patch assembly which contains some genes (e.g. TRBJ1) which are missing from the main assembly.
        // The issue does not exist in GRCh38 as that assembly is more complete.
        var refAlignments = runBwaMem(
            querySequences, mRefGenomeDictPath, mRefGenomeBwaIndexImagePath, ANNOTATION_ALIGN_SCORE_MIN, numThreads)
        if(mRefGenomeVersion == RefGenomeVersion.V37)
        {
            val patchAlignments = runGRCh37PatchAlignment(querySequences, ANNOTATION_ALIGN_SCORE_MIN, numThreads)
            require(patchAlignments.size == refAlignments.size)
            refAlignments = refAlignments.zip(patchAlignments).map { it.first + it.second }
        }

        val imgtAlignments = runImgtAlignment(mImgtSequences, querySequences, ANNOTATION_ALIGN_SCORE_MIN, numThreads)

        val annotations = processAlignments(alignmentMetadata, refAlignments, imgtAlignments)

        AlignmentMatchTsvWriter.write(outputDir, sampleId, annotations)

        return annotations
    }

    private data class AlignmentMetadata(val vdj: VDJSequence, val querySeqRange: IntRange, val querySeq: String)

    private fun processAlignments(
        metadatas: List<AlignmentMetadata>, refAlignments: List<List<Alignment>>, imgtAlignments: List<List<Alignment>>
    ): Collection<AlignmentAnnotation>
    {
        sLogger.debug("Processing alignments")
        val alignmentAnnotations = ArrayList<AlignmentAnnotation>()
        require(metadatas.size == refAlignments.size)
        require(metadatas.size == imgtAlignments.size)
        for ((index, metadata) in metadatas.withIndex())
        {
            alignmentAnnotations.add(processAlignments(metadata, refAlignments[index], imgtAlignments[index]))
        }
        sLogger.debug("Done processing alignments")
        return alignmentAnnotations
    }

    private fun processAlignments(
        metadata: AlignmentMetadata, refAlignments: Collection<Alignment>, imgtAlignments: Collection<Alignment>): AlignmentAnnotation
    {
        tryMatchReference(metadata, refAlignments)?.let { return it }
        return matchRearrangedGenes(metadata, imgtAlignments)
    }

    // Check if any alignments cover the whole sequence, in which case there is no VDJ rearrangement.
    private fun tryMatchReference(metadata: AlignmentMetadata, refAlignments: Collection<Alignment>): AlignmentAnnotation?
    {
        val alignments = preprocessAlignments(metadata, refAlignments)
        for (alignment in alignments)
        {
            val alignLength = (alignment.queryEnd - alignment.queryStart + 1)
            val minAlignLength = metadata.querySeq.length * ANNOTATION_MATCH_REF_IDENTITY
            if (alignment.editDistancePct <= 1 - ANNOTATION_MATCH_REF_IDENTITY && alignLength >= minAlignLength)
            {
                return AlignmentAnnotation(metadata.vdj, AlignmentStatus.NO_REARRANGEMENT, fullAlignment = alignment)
            }
        }
        return null
    }

    private data class GeneMatchCandidate(
        val alignment: Alignment,
        val imgtSequence: ImgtSequenceFile.Sequence,
        val gene: IgTcrGene
    )

    // Annotate the individual genes which have been rearranged to form this sequence.
    private fun matchRearrangedGenes(metadata: AlignmentMetadata, imgtAlignments: Collection<Alignment>): AlignmentAnnotation
    {
        val alignments = preprocessAlignments(metadata, imgtAlignments)

        val vdjSequence = metadata.vdj

        // we freeze the locus here. Reason is that there are cases where a low identity match (92%) from another
        // locus supercedes a 100% identity match from the correct locus
        val locus: IgTcrLocus? = vdjSequence.vAnchor?.geneType?.locus ?: vdjSequence.jAnchor?.geneType?.locus
        require(locus != null)

        // Find candidate gene matches, from which we will then select the best and supplementary matches.
        val candidates = HashMap<IgTcrRegion, ArrayList<GeneMatchCandidate>>()
        for (alignment in alignments)
        {
            val imgtSequence = mImgtSequences.sequencesByContig[alignment.refContig]!!
            val vdjGene = mVdjGenes[imgtSequence.geneAllele] ?: continue

            // For V/J gene segments, require 90% identity as a baseline.
            if (vdjGene.region.isVJ && alignment.editDistancePct > 1 - ANNOTATION_VJ_IDENTITY_MIN)
            {
                continue
            }

            // we must check to make sure the locus matches the top alignment
            // this ensure we do not annotate incorrect genes
            val geneLocus = IgTcrLocus.fromGeneName(vdjGene.geneName)
            if (locus != geneLocus)
            {
                continue
            }

            val candidate = GeneMatchCandidate(alignment, imgtSequence, vdjGene)
            candidates.getOrPut(vdjGene.region) { ArrayList() }.add(candidate)
        }

        val geneMatches = candidates.mapValues { selectBestGene(it.value) }
        val geneAnnotations = geneMatches.mapValues { (_, m) -> m?.let { m ->
            GeneAnnotation(
                m.primary.gene,
                m.primary.alignment,
                // For now, we are ignoring indels. May want to reevaluate that in the future.
                compareGeneSequence(metadata.querySeq, m.primary.alignment).alignedPctIdentity,
                m.supplementary.map { c -> c.gene })
        } }
        val alignmentStatus = getAlignmentStatus(geneMatches)

        return AlignmentAnnotation(
            vdjSequence = vdjSequence,
            status = alignmentStatus,
            vGene = geneAnnotations[IgTcrRegion.V_REGION],
            dGene = geneAnnotations[IgTcrRegion.D_REGION],
            jGene = geneAnnotations[IgTcrRegion.J_REGION]
        )
    }

    private fun compareSequenceToImgt(metadata: AlignmentMetadata, match: GeneMatchCandidate): AlignedSeqCompare
    {
        // Calculate percentage identity of the V region between the sample and IMGT sequence.
        // This is a heuristic for the degree of somatic hypermutation, which is a prognostic indicator for chronic lymphocytic leukemia.
        // https://pmc.ncbi.nlm.nih.gov/articles/PMC7248390/
        // The V region is from Cys104 (last amino acid of anchor) upstream to the start of the V exon.

        // TODO: want to compare the V side from Cys104 upstream until the IMGT sequence ends
        // TODO: want to compare the J side from anchor downstream until the IMGT sequence ends
        val refSeq = imgtSequence.sequenceWithRef.substring(alignment.refStart - 1, alignment.refEnd)
        val queryAlignedSeq = if (alignment.strand == Strand.FORWARD) querySeq else reverseComplement(querySeq)
        return compareAlignedSequence(queryAlignedSeq, refSeq, alignment.cigar)
    }

    private companion object
    {
        fun alignmentQuerySeqRange(vdj: VDJSequence) : IntRange
        {
            val fullSeq = vdj.layout.consensusSequenceString()
            val vdjSeq = vdj.sequence
            return alignmentQuerySeqRange(vdjSeq, fullSeq)
        }

        fun alignmentQuerySeqRange(vdjSeq: String, fullSeq: String): IntRange
        {
            val i = fullSeq.indexOf(vdjSeq)
            require(i >= 0)
            val range = max(i - ANNOTATION_VDJ_FLANK_BASES, 0) until
                    min(i + vdjSeq.length + ANNOTATION_VDJ_FLANK_BASES, fullSeq.length)
            return range
        }

        fun preprocessAlignments(metadata: AlignmentMetadata, alignments: Collection<Alignment>): List<Alignment>
        {
            val alignStartOffset = metadata.querySeqRange.start
            return alignments
                .filter { it.alignmentScore >= ANNOTATION_ALIGN_SCORE_MIN }
                .map {
                    // we also need to fix up the matches, since we did not use the full sequence to query, the queryAlignStart
                    // and queryAlignEnd indices are off
                    it.copy(
                        queryStart = it.queryStart + alignStartOffset,
                        queryEnd = it.queryEnd + alignStartOffset)
                }
                .sortedBy { -it.alignmentScore }
        }

        data class GeneMatch(
            val primary: GeneMatchCandidate,
            val supplementary: Collection<GeneMatchCandidate>
        )

        fun selectBestGene(candidateGenes: Collection<GeneMatchCandidate>) : GeneMatch?
        {
            if (candidateGenes.isEmpty())
            {
                return null
            }

            // Find the best gene match. There may be multiple, in which case we use a tie breaker and select the rest as supplementary matches.
            val baseComparator = Comparator(::compareGeneMatch)
            val tieBreakerComparator = Comparator(::geneMatchTieBreaker)
            val primary = candidateGenes.minWith(baseComparator.thenComparing(tieBreakerComparator))
            val supplementary = candidateGenes.filter { it != primary && baseComparator.compare(it, primary) == 0 }
            return GeneMatch(primary, supplementary)
        }

        // Better gene match compares less than.
        fun compareGeneMatch(candidate1: GeneMatchCandidate, candidate2: GeneMatchCandidate): Int
        {
            // Always prefer higher alignment score
            if (candidate1.alignment.alignmentScore > candidate2.alignment.alignmentScore)
            {
                return -1
            }
            else if (candidate1.alignment.alignmentScore < candidate2.alignment.alignmentScore)
            {
                return 1
            }

            // if scores are equal, we prefer the functional one
            if (candidate1.gene.isFunctional && !candidate2.gene.isFunctional)
            {
                return -1
            }
            else if (!candidate1.gene.isFunctional && candidate2.gene.isFunctional)
            {
                return 1
            }

            return 0
        }

        // Deterministic tie breaker for otherwise identical gene alignments
        fun geneMatchTieBreaker(candidate1: GeneMatchCandidate, candidate2: GeneMatchCandidate) : Int
        {
            return candidate1.gene.geneAllele.compareTo(candidate2.gene.geneAllele)
        }

        data class AlignedSeqCompare(
            val alignedBases: Int,
            val unalignedBases: Int,
            val matches: Int
        )
        {
            val alignedPctIdentity: Double = matches.toDouble() / alignedBases
        }

        fun compareAlignedSequence(query: String, ref: String, cigar: List<CigarElement>): AlignedSeqCompare
        {
            var alignedBases = 0
            var unalignedBases = 0
            var matches = 0

            var refIndex = 0
            var queryIndex = 0

            for (element in cigar)
            {
                val op = element.operator
                if (op.isAlignment)
                {
                    for (i in 0 until element.length)
                    {
                        if (query[queryIndex] == ref[refIndex])
                        {
                            matches++
                        }
                        queryIndex += 1
                        refIndex += 1
                    }
                    alignedBases += element.length
                }
                else if (op.consumesReadBases() || op.isClipping)
                {
                    unalignedBases += element.length
                    queryIndex += element.length
                }
                else if (op.consumesReferenceBases())
                {
                    unalignedBases += element.length
                    refIndex += element.length
                }
            }
            require(queryIndex == query.length)

            return AlignedSeqCompare(alignedBases, unalignedBases, matches)
        }

        fun getAlignmentStatus(geneMatches: Map<IgTcrRegion, GeneMatch?>): AlignmentStatus
        {
            val vMatched = geneMatches[IgTcrRegion.V_REGION] != null
            val dMatched = geneMatches[IgTcrRegion.D_REGION] != null
            val jMatched = geneMatches[IgTcrRegion.J_REGION] != null
            return if (vMatched)
                    (if (dMatched)
                            (if (jMatched) AlignmentStatus.V_D_J else AlignmentStatus.V_D)
                    else
                            (if (jMatched) AlignmentStatus.V_J else AlignmentStatus.V_ONLY))
            else
                    (if (dMatched)
                            (if (jMatched) AlignmentStatus.D_J else AlignmentStatus.D_ONLY)
                    else
                            (if (jMatched) AlignmentStatus.J_ONLY else AlignmentStatus.NO_VDJ_ALIGNMENT))
        }
    }
}
