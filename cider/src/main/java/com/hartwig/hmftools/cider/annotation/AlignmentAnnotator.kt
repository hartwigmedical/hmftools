package com.hartwig.hmftools.cider.annotation

import com.google.common.collect.Multimap
import com.google.common.collect.Multimaps
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.IgTcrGene.Companion.fromCommonIgTcrGene
import com.hartwig.hmftools.cider.AlignmentUtil.parseChromosome
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import org.apache.logging.log4j.LogManager
import java.util.*

enum class AlignmentStatus
{
    V_D_J, V_J, V_D, D_J, V_ONLY, D_ONLY, J_ONLY, NO_REARRANGEMENT, NO_VDJ_ALIGNMENT, SKIPPED_ALIGN
}

data class AlignmentAnnotation(
    val vdjSequence: VDJSequence,
    val vGene: IgTcrGene? = null,
    val vGeneSupplementary: List<IgTcrGene> = emptyList(),
    val dGene: IgTcrGene? = null,
    val dGeneSupplementary: List<IgTcrGene> = emptyList(),
    val jGene: IgTcrGene? = null,
    val jGeneSupplementary: List<IgTcrGene> = emptyList(),
    val vAlignment: AlignmentUtil.BwaMemAlignment? = null,
    val dAlignment: AlignmentUtil.BwaMemAlignment? = null,
    val jAlignment: AlignmentUtil.BwaMemAlignment? = null,
    val fullAlignment: AlignmentUtil.BwaMemAlignment? = null,
    val alignmentStatus: AlignmentStatus)

class AlignmentAnnotator
{
    private val sLogger = LogManager.getLogger(AlignmentAnnotator::class.java)

    private val mRefGenome: RefGenomeSource
    private val mRefGenomeBwaIndexImagePath: String
    private val mVdjGenes: Map<Pair<String, Strand>, List<IgTcrGene>>

    // class to help associate the data back
    data class AlignmentRunData(val vdj: VDJSequence, val key: Int, val querySeqRange: IntRange, val querySeq: String)

    constructor(refGenomeFastaPath: String, refGenomeBwaIndexImagePath: String)
    {
        this.mRefGenome = RefGenomeSource.loadRefGenome(refGenomeFastaPath)
        this.mRefGenomeBwaIndexImagePath = refGenomeBwaIndexImagePath

        val refGenomeVersion = deriveRefGenomeVersion(mRefGenome)

        // Explicitly using an ArrayList here to give a deterministic iteration order when finding genes later.
        val vdjGenes: HashMap<Pair<String, Strand>, ArrayList<IgTcrGene>> = HashMap()

        val igTcrGenes = IgTcrGeneFile.read(refGenomeVersion)
            .map { o -> fromCommonIgTcrGene(o) }

        // find all the genes that we need
        for (geneData in igTcrGenes)
        {
            if (geneData.region !in arrayOf(IgTcrRegion.V_REGION, IgTcrRegion.D_REGION, IgTcrRegion.J_REGION))
            {
                continue
            }

            if (geneData.geneLocation == null)
            {
                continue
            }

            val key = Pair(geneData.geneLocation.chromosome, geneData.geneLocation.strand)
            vdjGenes.computeIfAbsent(key) { ArrayList() }.add(geneData)

            /* sLogger.debug(
                "found constant region gene: {}, type: {}, location: {}",
                geneData.GeneName, igConstantRegionType, genomeRegionStrand
            )*/
        }

        mVdjGenes = vdjGenes
    }

    fun runAnnotate(sampleId: String, vdjList: List<VDJSequence>, outputDir: String, numThreads: Int)
            : Collection<AlignmentAnnotation>
    {
        sLogger.info("Running alignment annotation")

        // assign a key to each VDJ, such that we can keep track of them
        var key = 0
        val alignmentRunDataMap : MutableMap<Int, AlignmentRunData> = HashMap()

        for (vdj in vdjList)
        {
            val querySeqRange = alignmentQuerySeqRange(vdj)
            val alignmentRunData = AlignmentRunData(vdj,
                key++,
                querySeqRange,
                vdj.layout.consensusSequenceString().substring(querySeqRange))
            alignmentRunDataMap[alignmentRunData.key] = alignmentRunData
        }

        val alignmentResults = AlignmentUtil.runBwaMem(
            alignmentRunDataMap.mapValues { runData -> runData.value.querySeq },
             mRefGenome, mRefGenomeBwaIndexImagePath, BWA_ALIGNMENT_SCORE_MIN, numThreads)

        // put all into an identity hash multimap
        val vdjToAlignment: Multimap<AlignmentRunData, AlignmentUtil.BwaMemAlignment> = Multimaps.newListMultimap(IdentityHashMap()) { ArrayList() }

        for ((vdjKey, alignment) in alignmentResults.entries())
        {
            val alignmentRunData = alignmentRunDataMap[vdjKey]

            if (alignmentRunData == null)
            {
                sLogger.fatal("error processing alignment results: cannot find key: {}", vdjKey)
                throw RuntimeException("error processing alignment results: cannot find key: $vdjKey")
            }

            vdjToAlignment.put(alignmentRunData, alignment)
        }

        val annotations = processAlignments(alignmentRunDataMap.values, vdjToAlignment)

        AlignmentMatchTsvWriter.write(outputDir, sampleId, annotations)

        return annotations
    }

    // process the alignment matches for each VDJ, and set the alignmentAnnotation in the VdjAnnotation
    // NOTE: we cannot use matches.keySet(), as it might not include some VDJs that returned no match
    fun processAlignments(alignmentRunDataList: Collection<AlignmentRunData>, alignments: Multimap<AlignmentRunData, AlignmentUtil.BwaMemAlignment>)
    : Collection<AlignmentAnnotation>
    {
        val alignmentAnnotations = ArrayList<AlignmentAnnotation>()
        for (runData in alignmentRunDataList)
        {
            alignmentAnnotations.add(processAlignments(runData, alignments[runData]))
        }
        return alignmentAnnotations
    }

    fun processAlignments(alignmentRunData: AlignmentRunData, alignments: Collection<AlignmentUtil.BwaMemAlignment>)
    : AlignmentAnnotation
    {
        val vdjSequence: VDJSequence = alignmentRunData.vdj

        val alignStartOffset = alignmentRunData.querySeqRange.start

        val alignments = alignments
            .filter { m -> m.alignmentScore >= BWA_ALIGNMENT_SCORE_MIN }
            .map {
                // we also need to fix up the matches, since we did not use the full sequence to query, the queryAlignStart
                // and queryAlignEnd indices are off
                it.copy(
                    queryAlignStart = it.queryAlignStart + alignStartOffset,
                    queryAlignEnd = it.queryAlignEnd + alignStartOffset)
            }
            .sortedBy { m -> -m.alignmentScore }

        // Check if any alignments cover the whole sequence, in which case there is no VDJ rearrangement.
        for (alignment in alignments)
        {
            if (alignment.percentageIdent >= CiderConstants.ALIGNMENT_MATCH_FULL_MATCH_IDENTITY &&
                alignmentRunData.querySeq.length <= (alignment.queryAlignEnd - alignment.queryAlignStart) + 5
            )
            {
                return AlignmentAnnotation(
                    vdjSequence = vdjSequence,
                    fullAlignment = alignment,
                    alignmentStatus = AlignmentStatus.NO_REARRANGEMENT
                )
            }
        }

        // we freeze the locus here. Reason is that there are cases where a low identity match (92%) from another
        // locus supercedes a 100% identity match from the correct locus
        val locus: IgTcrLocus? = vdjSequence.vAnchor?.geneType?.locus ?: vdjSequence.jAnchor?.geneType?.locus
        require(locus != null)

        // Find candidate gene matches, from which we will then select the best and supplementary matches.
        val vGeneCandidates: MutableList<Pair<IgTcrGene, AlignmentUtil.BwaMemAlignment>> = ArrayList()
        val dGeneCandidates: MutableList<Pair<IgTcrGene, AlignmentUtil.BwaMemAlignment>> = ArrayList()
        val jGeneCandidates: MutableList<Pair<IgTcrGene, AlignmentUtil.BwaMemAlignment>> = ArrayList()
        for (alignment in alignments)
        {
            val vdjGene: IgTcrGene? = findGene(alignment)
            if (vdjGene == null)
            {
                continue
            }

            // for V/J gene segments, we mandate 90% identity
            if (vdjGene.region in arrayOf(IgTcrRegion.V_REGION, IgTcrRegion.J_REGION) && alignment.percentageIdent < CiderConstants.ALIGNMENT_MATCH_MIN_VJ_IDENTITY)
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
            vGene = vGeneMatch?.gene,
            dGene = dGeneMatch?.gene,
            jGene = jGeneMatch?.gene,
            vAlignment = vGeneMatch?.alignment,
            dAlignment = dGeneMatch?.alignment,
            jAlignment = jGeneMatch?.alignment,
            vGeneSupplementary = vGeneMatch?.supplementaryGenes ?: emptyList(),
            dGeneSupplementary = dGeneMatch?.supplementaryGenes ?: emptyList(),
            jGeneSupplementary = jGeneMatch?.supplementaryGenes ?: emptyList(),
            alignmentStatus = alignmentStatus)
    }

    fun findGene(alignment: AlignmentUtil.BwaMemAlignment) : IgTcrGene?
    {
        val chromosome = parseChromosome(alignment.refContig)

        val geneDataList = mVdjGenes[Pair(chromosome, alignment.strand)] ?: return null

        var bestGene : IgTcrGene? = null

        for (gene in geneDataList)
        {
            val geneLocation = gene.geneLocation ?: continue

            require(geneLocation.chromosome == chromosome)
            require(geneLocation.strand == alignment.strand)

            if (bestGene == null || !bestGene.isFunctional)
            {
                // check if they overlap. We prioritise functional genes
                if (geneLocation.posStart <= alignment.refEnd + GENE_REGION_TOLERANCE &&
                    geneLocation.posEnd >= alignment.refStart - GENE_REGION_TOLERANCE)
                {
                    bestGene = gene
                }
            }
        }

        return bestGene
    }

    companion object
    {
        // Require a match of minimum ~20 bases. If we want to match D segment that is shorter
        // we will need a higher cut off, maybe 10, but will get many false positive hits that are longer but more mismatches
        const val BWA_ALIGNMENT_SCORE_MIN = 20

        const val FLANKING_BASES = 50

        // D genes often are very short, for example, TRBD1 is only 12 bases. We allow more leeway to match
        // an alignment to the gene
        const val GENE_REGION_TOLERANCE = 50

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
            val range = Math.max(i - FLANKING_BASES, 0) until
                    Math.min(i + vdjSeq.length + FLANKING_BASES, fullSeq.length)
            return range
        }

        data class GeneMatch(val gene: IgTcrGene, val alignment: AlignmentUtil.BwaMemAlignment, val supplementaryGenes: List<IgTcrGene>)

        fun selectBestGene(candidateGenes: List<Pair<IgTcrGene, AlignmentUtil.BwaMemAlignment>>) : GeneMatch?
        {
            if (candidateGenes.isEmpty())
            {
                return null
            }

            // Find the best gene match. There may be multiple, in which case we use a tie breaker and select the rest as supplementary matches.
            val baseComparator = Comparator<Pair<IgTcrGene, AlignmentUtil.BwaMemAlignment>>
                { p1, p2 -> compareGeneMatch(p1.second, p1.first, p2.second, p2.first) }
            val tieBreakerComparator = Comparator<Pair<IgTcrGene, AlignmentUtil.BwaMemAlignment>>
                { p1, p2 -> geneMatchTieBreaker(p1.first, p2.first) }
            val best = candidateGenes.minWith(baseComparator.thenComparing(tieBreakerComparator))
            val supplementary = candidateGenes.filter { it != best && baseComparator.compare(it, best) == 0 }.map { it.first }
            return GeneMatch(best.first, best.second, supplementary)
        }

        // Better gene match compares less than.
        fun compareGeneMatch(alignment1: AlignmentUtil.BwaMemAlignment, gene1: IgTcrGene, alignment2: AlignmentUtil.BwaMemAlignment, gene2: IgTcrGene)
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
            return gene1.geneName.compareTo(gene2.geneName)
        }
    }
}
