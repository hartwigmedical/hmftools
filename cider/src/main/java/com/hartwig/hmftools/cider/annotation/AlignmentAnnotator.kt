package com.hartwig.hmftools.cider.annotation

import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_MATCH_REF_IDENTITY
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_ALIGN_SCORE_MIN
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_VDJ_FLANK_BASES
import com.hartwig.hmftools.cider.CiderConstants.ANNOTATION_VJ_IDENTITY_MIN
import com.hartwig.hmftools.cider.CiderUtils.getResourceAsFile
import com.hartwig.hmftools.cider.CiderUtils.getResourceAsStream
import com.hartwig.hmftools.cider.genes.IgTcrGene
import com.hartwig.hmftools.cider.genes.IgTcrGene.Companion.fromCommonIgTcrGene
import com.hartwig.hmftools.cider.genes.IgTcrLocus
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
import com.hartwig.hmftools.common.cider.IgTcrRegion
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
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
    val comparison: ShmGeneComparison?,
    val supplementaryGenes: List<IgTcrGene>
)

data class AlignmentAnnotation(
    val vdjSequence: VDJSequence,
    val alignmentQueryRange: IntRange,
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
        val alignments = preprocessAlignments(refAlignments)
        for (alignment in alignments)
        {
            val distance = alignment.editDistance.toDouble() / metadata.querySeq.length
            val minAlignLength = metadata.querySeq.length * ANNOTATION_MATCH_REF_IDENTITY
            if (distance <= 1 - ANNOTATION_MATCH_REF_IDENTITY && alignment.queryAlignLength >= minAlignLength)
            {
                return AlignmentAnnotation(
                    metadata.vdj,
                    metadata.querySeqRange,
                    AlignmentStatus.NO_REARRANGEMENT,
                    fullAlignment = alignment)
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
        val alignments = preprocessAlignments(imgtAlignments)

        val vdj = metadata.vdj

        // we freeze the locus here. Reason is that there are cases where a low identity match (92%) from another
        // locus supercedes a 100% identity match from the correct locus
        val locus: IgTcrLocus? = vdj.vAnchor?.geneType?.locus ?: vdj.jAnchor?.geneType?.locus
        require(locus != null)

        // Find candidate gene matches, from which we will then select the best and supplementary matches.
        val candidates = HashMap<IgTcrRegion, ArrayList<GeneMatchCandidate>>()
        for (alignment in alignments)
        {
            val imgtSequence = mImgtSequences.sequencesByContig[alignment.refContig]!!
            val vdjGene = mVdjGenes[imgtSequence.geneAllele] ?: continue

            if (vdjGene.region.isVJ)
            {
                // For V/J gene segments, require 90% identity as a baseline.
                val distance = alignment.alignedEditDistance.toDouble() / alignment.queryAlignLength
                if (distance > 1 - ANNOTATION_VJ_IDENTITY_MIN)
                {
                    continue
                }

                // Also require that the alignment covers the anchor.
                // TODO: what if the anchor matches the ref genome part?
                val layoutAlignRange = alignment.queryRange.start + metadata.querySeqRange.start..alignment.queryRange.endInclusive + metadata.querySeqRange.start
                val anchorBoundary = if (vdjGene.region == IgTcrRegion.V_REGION) vdj.vAnchorBoundary?.plus(vdj.layoutSliceStart - 1)
                    else vdj.jAnchorBoundary?.plus(vdj.layoutSliceStart)
                if (anchorBoundary != null && !layoutAlignRange.contains(anchorBoundary))
                {
                    continue
                }
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
                compareSequenceToImgt(metadata, m.primary),
                m.supplementary.map { c -> c.gene })
        } }
        val alignmentStatus = getAlignmentStatus(geneMatches)

        return AlignmentAnnotation(
            vdjSequence = vdj,
            alignmentQueryRange = metadata.querySeqRange,
            status = alignmentStatus,
            vGene = geneAnnotations[IgTcrRegion.V_REGION],
            dGene = geneAnnotations[IgTcrRegion.D_REGION],
            jGene = geneAnnotations[IgTcrRegion.J_REGION]
        )
    }

    private fun compareSequenceToImgt(metadata: AlignmentMetadata, match: GeneMatchCandidate): ShmGeneComparison?
    {
        return when(match.gene.region)
        {
            IgTcrRegion.V_REGION -> compareVRegionToImgt(metadata.vdj, match.imgtSequence, metadata.querySeqRange, match.alignment)
            IgTcrRegion.J_REGION -> compareJRegionToImgt(metadata.vdj, match.imgtSequence, metadata.querySeqRange, match.alignment)
            else -> null
        }
    }

    private companion object
    {
        private val sLogger = LogManager.getLogger(AlignmentAnnotator::class.java)

        fun alignmentQuerySeqRange(vdj: VDJSequence) : IntRange
        {
            val fullSeq = vdj.layout.consensusSequenceString()
            val vdjSeq = vdj.sequence
            val i = fullSeq.indexOf(vdjSeq)
            require(i >= 0)
            val range = max(i - ANNOTATION_VDJ_FLANK_BASES, 0) until
                    min(i + vdjSeq.length + ANNOTATION_VDJ_FLANK_BASES, fullSeq.length)
            return range
        }

        // Run alignment against a patch of the GRCh37 genome which includes more genes, particularly TRBJ1.
        fun runGRCh37PatchAlignment(sequences: List<String>, alignScoreThreshold: Int, threadCount: Int): List<List<Alignment>>
        {
            sLogger.debug("Aligning ${sequences.size} sequences to GRCh37 patch")
            val refFastaName = "7_gl582971_fix.fasta"
            val refDictStream = getResourceAsStream("$refFastaName.dict")
            val refIndexImagePath = getResourceAsFile("$refFastaName.img")
            return runBwaMem(sequences, refDictStream, refIndexImagePath, alignScoreThreshold, threadCount)
        }

        // Run alignment against the IMGT gene allele sequences.
        fun runImgtAlignment(imgtSequenceFile: ImgtSequenceFile, sequences: List<String>, alignScoreThreshold: Int, threadCount: Int): List<List<Alignment>>
        {
            sLogger.debug("Aligning ${sequences.size} sequences to IMGT gene alleles")
            return runBwaMem(sequences, imgtSequenceFile.fastaDict, imgtSequenceFile.bwamemImgPath, alignScoreThreshold, threadCount)
        }

        fun preprocessAlignments(alignments: Collection<Alignment>): List<Alignment>
        {
            return alignments
                .filter { it.alignmentScore >= ANNOTATION_ALIGN_SCORE_MIN }
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
            // Always prefer higher alignment score.
            if (candidate1.alignment.alignmentScore > candidate2.alignment.alignmentScore)
            {
                return -1
            }
            else if (candidate1.alignment.alignmentScore < candidate2.alignment.alignmentScore)
            {
                return 1
            }

            return 0
        }

        // Deterministic tie breaker for otherwise identical gene alignments
        fun geneMatchTieBreaker(candidate1: GeneMatchCandidate, candidate2: GeneMatchCandidate) : Int
        {
            // Prefer the functional gene.
            if (candidate1.gene.isFunctional && !candidate2.gene.isFunctional)
            {
                return -1
            }
            else if (!candidate1.gene.isFunctional && candidate2.gene.isFunctional)
            {
                return 1
            }

            // Use gene name as a final tie-breaker.
            return candidate1.gene.geneAllele.compareTo(candidate2.gene.geneAllele)
        }

        fun compareVRegionToImgt(
            vdj: VDJSequence, imgtSequence: ImgtSequenceFile.Sequence, queryRange: IntRange, alignment: Alignment) : ShmGeneComparison?
        {
            val layoutAnchorBoundary = vdj.layoutSliceStart + (vdj.vAnchorBoundary ?: return null)
            return compareVJRegionToImgt(
                vdj.layout.consensusSequenceString(), VJ.V, layoutAnchorBoundary, imgtSequence, queryRange, alignment)
        }

        fun compareJRegionToImgt(
            vdj: VDJSequence, imgtSequence: ImgtSequenceFile.Sequence, queryRange: IntRange, alignment: Alignment) : ShmGeneComparison?
        {
            val layoutAnchorBoundary = vdj.layoutSliceStart + (vdj.jAnchorBoundary ?: return null)
            return compareVJRegionToImgt(
                vdj.layout.consensusSequenceString(), VJ.J, layoutAnchorBoundary, imgtSequence, queryRange, alignment)
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
