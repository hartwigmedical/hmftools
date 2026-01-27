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
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.math.max
import kotlin.math.min

enum class AlignmentStatus
{
    V_D_J, V_J, V_D, D_J, V_ONLY, D_ONLY, J_ONLY, NO_REARRANGEMENT, NO_VDJ_ALIGNMENT, SKIPPED_ALIGN
}

data class AlignmentAnnotation(
    val vdjSequence: VDJSequence,
    val status: AlignmentStatus,
    val fullAlignment: Alignment? = null,
    val vGene: IgTcrGene? = null,
    val vAlignment: Alignment? = null,
    val vIdentity: Double? = null,
    val vGeneSupplementary: List<IgTcrGene> = emptyList(),
    val dGene: IgTcrGene? = null,
    val dAlignment: Alignment? = null,
    val dIdentity: Double? = null,
    val dGeneSupplementary: List<IgTcrGene> = emptyList(),
    val jGene: IgTcrGene? = null,
    val jAlignment: Alignment? = null,
    val jIdentity: Double? = null,
    val jGeneSupplementary: List<IgTcrGene> = emptyList()
)

class AlignmentAnnotator
{
    private val sLogger = LogManager.getLogger(AlignmentAnnotator::class.java)

    private val mRefGenomeVersion: RefGenomeVersion
    private val mRefGenomeDictPath: String
    private val mRefGenomeBwaIndexImagePath: String
    private val mVdjGenes: Map<String, IgTcrGene>
    private val mImgtSequences: ImgtSequenceFile

    private data class AlignmentMetadata(val vdj: VDJSequence, val querySeqRange: IntRange, val querySeq: String)

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
        val vGeneCandidates = ArrayList<Pair<IgTcrGene, Alignment>>()
        val dGeneCandidates = ArrayList<Pair<IgTcrGene, Alignment>>()
        val jGeneCandidates = ArrayList<Pair<IgTcrGene, Alignment>>()
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

            if (vdjGene.region == IgTcrRegion.V_REGION)
            {
                vGeneCandidates.add(Pair(vdjGene, alignment))
            }
            if (vdjGene.region == IgTcrRegion.D_REGION)
            {
                dGeneCandidates.add(Pair(vdjGene, alignment))
            }
            if (vdjGene.region == IgTcrRegion.J_REGION)
            {
                jGeneCandidates.add(Pair(vdjGene, alignment))
            }
        }

        val vGeneMatch = selectBestGene(vGeneCandidates)
        val dGeneMatch = selectBestGene(dGeneCandidates)
        val jGeneMatch = selectBestGene(jGeneCandidates)

        // determine status
        val alignmentStatus: AlignmentStatus = if (vGeneMatch != null)
        {
            if (dGeneMatch != null)
            {
                if (jGeneMatch != null)
                {
                    AlignmentStatus.V_D_J
                }
                else
                {
                    AlignmentStatus.V_D
                }
            }
            else if (jGeneMatch != null)
            {
                AlignmentStatus.V_J
            }
            else
            {
                AlignmentStatus.V_ONLY
            }
        }
        else
        {
            if (dGeneMatch != null)
            {
                if (jGeneMatch != null)
                {
                    AlignmentStatus.D_J
                }
                else
                {
                    AlignmentStatus.D_ONLY
                }
            }
            else if (jGeneMatch != null)
            {
                AlignmentStatus.J_ONLY
            }
            else
            {
                AlignmentStatus.NO_VDJ_ALIGNMENT
            }
        }

        return AlignmentAnnotation(
            vdjSequence = vdjSequence,
            status = alignmentStatus,
            vGene = vGeneMatch?.gene,
            dGene = dGeneMatch?.gene,
            jGene = jGeneMatch?.gene,
            vAlignment = vGeneMatch?.alignment,
            dAlignment = dGeneMatch?.alignment,
            jAlignment = jGeneMatch?.alignment,
            vGeneSupplementary = vGeneMatch?.supplementaryGenes ?: emptyList(),
            dGeneSupplementary = dGeneMatch?.supplementaryGenes ?: emptyList(),
            jGeneSupplementary = jGeneMatch?.supplementaryGenes ?: emptyList())
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

        data class GeneMatch(val gene: IgTcrGene, val alignment: Alignment, val supplementaryGenes: List<IgTcrGene>)

        fun selectBestGene(candidateGenes: List<Pair<IgTcrGene, Alignment>>) : GeneMatch?
        {
            if (candidateGenes.isEmpty())
            {
                return null
            }

            // Find the best gene match. There may be multiple, in which case we use a tie breaker and select the rest as supplementary matches.
            val baseComparator = Comparator<Pair<IgTcrGene, Alignment>>
                { p1, p2 -> compareGeneMatch(p1.second, p1.first, p2.second, p2.first) }
            val tieBreakerComparator = Comparator<Pair<IgTcrGene, Alignment>>
                { p1, p2 -> geneMatchTieBreaker(p1.first, p2.first) }
            val best = candidateGenes.minWith(baseComparator.thenComparing(tieBreakerComparator))
            val supplementary = candidateGenes.filter { it != best && baseComparator.compare(it, best) == 0 }.map { it.first }
            return GeneMatch(best.first, best.second, supplementary)
        }

        // Better gene match compares less than.
        fun compareGeneMatch(alignment1: Alignment, gene1: IgTcrGene, alignment2: Alignment, gene2: IgTcrGene)
            : Int
        {
            // Always prefer higher alignment score
            if (alignment1.alignmentScore > alignment2.alignmentScore)
            {
                return -1
            }
            else if (alignment1.alignmentScore < alignment2.alignmentScore)
            {
                return 1
            }

            // if scores are equal, we prefer the functional one
            if (gene1.isFunctional && !gene2.isFunctional)
            {
                return -1
            }
            else if (!gene1.isFunctional && gene2.isFunctional)
            {
                return 1
            }

            return 0
        }

        // Deterministic tie breaker for otherwise identical gene alignments
        fun geneMatchTieBreaker(gene1: IgTcrGene, gene2: IgTcrGene) : Int
        {
            return gene1.geneAllele.compareTo(gene2.geneAllele)
        }
    }
}
