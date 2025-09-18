package com.hartwig.hmftools.cider.alignment

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.google.common.collect.Multimaps
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.IgTcrGene.Companion.fromCommonIgTcrGene
import com.hartwig.hmftools.cider.alignment.AlignmentUtil.parseChromosome
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
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
    val dGene: IgTcrGene? = null,
    val jGene: IgTcrGene? = null,
    val vMatch: AlignmentUtil.BwaMemMatch? = null,
    val dMatch: AlignmentUtil.BwaMemMatch? = null,
    val jMatch: AlignmentUtil.BwaMemMatch? = null,
    val fullMatch: AlignmentUtil.BwaMemMatch? = null,
    val alignmentStatus: AlignmentStatus)

class AlignmentAnnotator
{
    private val sLogger = LogManager.getLogger(AlignmentAnnotator::class.java)

    val vdjGenes: Multimap<Pair<String, Strand>, IgTcrGene>

    // class to help associate the data back
    data class AlignmentRunData(val vdj: VDJSequence, val key: Int, val querySeqRange: IntRange, val querySeq: String)

    init
    {
        vdjGenes = ArrayListMultimap.create()

        // since we use V38 for alignment annotation
        val igTcrGenes = IgTcrGeneFile.read(CiderConstants.BLAST_REF_GENOME_VERSION)
            .map { o -> fromCommonIgTcrGene(o) }

        // find all the genes that are we need
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

            vdjGenes.put(Pair(geneData.geneLocation.chromosome, geneData.geneLocation.strand), geneData)

            /* sLogger.debug(
                "found constant region gene: {}, type: {}, location: {}",
                geneData.GeneName, igConstantRegionType, genomeRegionStrand
            )*/
        }
    }

    fun runAnnotate(sampleId: String, vdjList: List<VDJSequence>, refGenomeFastaPath: String, refGenomeIndexPath: String, outputDir: String, numThreads: Int)
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
             refGenomeFastaPath, refGenomeIndexPath, numThreads)

        // put all into an identity hash multimap
        val vdjToMatch: Multimap<AlignmentRunData, AlignmentUtil.BwaMemMatch> = Multimaps.newListMultimap(IdentityHashMap()) { ArrayList() }

        for ((vdjKey, match) in alignmentResults.entries())
        {
            val alignmentRunData = alignmentRunDataMap[vdjKey]

            if (alignmentRunData == null)
            {
                sLogger.fatal("error processing alignment results: cannot find key: {}", vdjKey)
                throw RuntimeException("error processing alignment results: cannot find key: $vdjKey")
            }

            vdjToMatch.put(alignmentRunData, match)
        }

        val annotations = processAlignmentMatches(alignmentRunDataMap.values, vdjToMatch)

        AlignmentMatchTsvWriter.write(outputDir, sampleId, annotations)

        return annotations
    }

    // process the alignment matches for each VDJ, and set the alignmentAnnotation in the VdjAnnotation
    // NOTE: we cannot use matches.keySet(), as it might not include some VDJs that returned no match
    fun processAlignmentMatches(alignmentRunDataList: Collection<AlignmentRunData>, matches: Multimap<AlignmentRunData, AlignmentUtil.BwaMemMatch>)
    : Collection<AlignmentAnnotation>
    {
        val alignmentAnnotations = ArrayList<AlignmentAnnotation>()
        for (runData in alignmentRunDataList)
        {
            alignmentAnnotations.add(processAlignmentMatches(runData, matches[runData]))
        }
        return alignmentAnnotations
    }

    fun processAlignmentMatches(alignmentRunData: AlignmentRunData, matches: Collection<AlignmentUtil.BwaMemMatch>)
    : AlignmentAnnotation
    {
        val vdjSequence: VDJSequence = alignmentRunData.vdj

        val alignStartOffset = alignmentRunData.querySeqRange.start

        val matches = matches
            .filter { m -> m.alignmentScore >= BWA_ALIGNMENT_SCORE_MIN }
            .map {
                // we also need to fix up the matches, since we did not use the full sequence to query, the queryAlignStart
                // and queryAlignEnd indices are off
                it.copy(
                    queryAlignStart = it.queryAlignStart + alignStartOffset,
                    queryAlignEnd = it.queryAlignEnd + alignStartOffset)
            }
            .sortedBy { m -> -m.alignmentScore }

        // we freeze the locus here. Reason is that there are cases where a low identity match (92%) from another
        // locus supercedes a 100% identity match from the correct locus
        val locus: IgTcrLocus? = vdjSequence.vAnchor?.geneType?.locus ?: vdjSequence.jAnchor?.geneType?.locus

        require(locus != null)

        var vGene: IgTcrGene? = null
        var vMatch: AlignmentUtil.BwaMemMatch? = null
        var dGene: IgTcrGene? = null
        var dMatch: AlignmentUtil.BwaMemMatch? = null
        var jGene: IgTcrGene? = null
        var jMatch: AlignmentUtil.BwaMemMatch? = null

        for (match in matches)
        {
            // check if it matches whole way. Require 95% sequence identity
            if (match.percentageIdent >= CiderConstants.ALIGNMENT_MATCH_FULL_MATCH_IDENTITY &&
                alignmentRunData.querySeq.length <= (match.queryAlignEnd - match.queryAlignStart) + 5)
            {
                // sLogger.debug("matches ref genome: {}", match.subjectTitle)
                // sLogger.debug("  query seq: {}", match.alignedPartOfQuerySeq)
                // sLogger.debug("subject seq: {}", match.alignedPartOfSubjectSeq)

                return AlignmentAnnotation(
                    vdjSequence = vdjSequence,
                    fullMatch = match,
                    alignmentStatus = AlignmentStatus.NO_REARRANGEMENT)
            }

            var vdjGene: IgTcrGene? = findGene(match)

            // for V/J gene segments, we mandate 90% identity
            if (vdjGene != null && vdjGene.region in arrayOf(IgTcrRegion.V_REGION, IgTcrRegion.J_REGION) && match.percentageIdent < CiderConstants.ALIGNMENT_MATCH_MIN_VJ_IDENTITY)
            {
                vdjGene = null
            }

            if (vdjGene != null)
            {
                val geneLocus = IgTcrLocus.fromGeneName(vdjGene.geneName)

                if (locus == geneLocus)
                {
                    // we must check to make sure the locus matches the top alignment
                    // this ensure we do not annotate incorrect genes
                    // we found a VDJ gene, see which one it is
                    // but also need to check against identity
                    if (vdjGene.region == IgTcrRegion.V_REGION && isBetterMatch(vMatch, vGene, match, vdjGene))
                    {
                        vGene = vdjGene
                        vMatch = match
                    }
                    if (vdjGene.region == IgTcrRegion.D_REGION && isBetterMatch(dMatch, dGene, match, vdjGene))
                    {
                        dGene = vdjGene
                        dMatch = match
                    }
                    if (vdjGene.region == IgTcrRegion.J_REGION && isBetterMatch(jMatch, jGene, match, vdjGene))
                    {
                        jGene = vdjGene
                        jMatch = match
                    }
                }
            }
        }

        // determine status
        val alignmentStatus: AlignmentStatus = if (vGene != null)
        {
            if (dGene != null)
            {
                if (jGene != null)
                {
                    AlignmentStatus.V_D_J
                }
                else
                {
                    AlignmentStatus.V_D
                }
            }
            else if (jGene != null)
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
            if (dGene != null)
            {
                if (jGene != null)
                {
                    AlignmentStatus.D_J
                }
                else
                {
                    AlignmentStatus.D_ONLY
                }
            }
            else if (jGene != null)
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
            vGene = vGene,
            dGene = dGene,
            jGene = jGene,
            vMatch = vMatch,
            dMatch = dMatch,
            jMatch = jMatch,
            alignmentStatus = alignmentStatus)
    }

    fun findGene(alignmentMatch: AlignmentUtil.BwaMemMatch) : IgTcrGene?
    {
        val chromosome = parseChromosome(alignmentMatch.refContig)

        val geneDataList = vdjGenes[Pair(chromosome, alignmentMatch.strand)]

        var bestGene : IgTcrGene? = null

        for (gene in geneDataList)
        {
            val geneLocation = gene.geneLocation ?: continue

            require(geneLocation.chromosome == chromosome)
            require(geneLocation.strand == alignmentMatch.strand)

            if (bestGene == null || !bestGene.isFunctional)
            {
                // check if they overlap. We prioritise functional genes
                if (geneLocation.posStart <= alignmentMatch.refEnd + GENE_REGION_TOLERANCE &&
                    geneLocation.posEnd >= alignmentMatch.refStart - GENE_REGION_TOLERANCE)
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

        private val sLogger = LogManager.getLogger(AlignmentAnnotator::class.java)

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

        fun isBetterMatch(existingMatch: AlignmentUtil.BwaMemMatch?, existingGene: IgTcrGene?, newMatch: AlignmentUtil.BwaMemMatch, newGene: IgTcrGene) : Boolean
        {
            if (existingMatch == null || existingGene == null)
            {
                return true
            }

            if (newMatch.alignmentScore > existingMatch.alignmentScore)
            {
                // always prefer higher alignment score
                return true
            }

            if ((newMatch.alignmentScore == existingMatch.alignmentScore))
            {
                if (!existingGene.isFunctional && newGene.isFunctional) {
                    // if scores are equal, we prefer the functional one
                    sLogger.trace(
                        "prefer functional gene: {}, alignScore: {} over non functional: {}, alignScore: {}",
                        newGene.geneAllele, newMatch.alignmentScore, existingGene.geneAllele, existingMatch.alignmentScore
                    )
                    return true
                }
                // Deterministic tie breaker for otherwise identical gene alignments
                if (newGene.geneName < existingGene.geneName) {
                    return true
                }
            }

            return false
        }
    }
}
