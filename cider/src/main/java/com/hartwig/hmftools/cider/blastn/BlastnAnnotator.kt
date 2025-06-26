package com.hartwig.hmftools.cider.blastn

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.google.common.collect.Multimaps
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.IgTcrGene.Companion.fromCommonIgTcrGene
import com.hartwig.hmftools.cider.blastn.BlastnUtil.toGenomicLocation
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.Doubles
import org.apache.logging.log4j.LogManager
import java.util.*

enum class BlastnStatus
{
    V_D_J, V_J, V_D, D_J, V_ONLY, D_ONLY, J_ONLY, NO_REARRANGEMENT, NO_VDJ_ALIGNMENT, SKIPPED_BLASTN
}

data class BlastnAnnotation(
    val vdjSequence: VDJSequence,
    val vGene: IgTcrGene? = null,
    val dGene: IgTcrGene? = null,
    val jGene: IgTcrGene? = null,
    val vMatch: BlastnMatch? = null,
    val dMatch: BlastnMatch? = null,
    val jMatch: BlastnMatch? = null,
    val fullMatch: BlastnMatch? = null,
    val blastnStatus: BlastnStatus)

// run blastn and matches and turn them into annotations
class BlastnAnnotator
{
    private val sLogger = LogManager.getLogger(BlastnAnnotator::class.java)

    val vdjGenes: Multimap<Pair<String, Strand>, IgTcrGene>

    // class to help associate the data back
    data class BlastnRunData(val vdj: VDJSequence, val key: Int, val querySeqRange: IntRange, val querySeq: String)

    init
    {
        vdjGenes = ArrayListMultimap.create()

        // since we use V38 for blast
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

    fun runAnnotate(sampleId: String, blastDir: String, blastDb: String, vdjList: List<VDJSequence>, outputDir: String, numThreads: Int)
            : Collection<BlastnAnnotation>
    {
        // assign a key to each VDJ, such that we can keep track of them
        var key = 0
        val blastnRunDataMap : MutableMap<Int, BlastnRunData> = HashMap()

        for (vdj in vdjList)
        {
            val querySeqRange = blastnQuerySeqRange(vdj)
            val blastnRunData = BlastnRunData(vdj,
                key++,
                querySeqRange,
                vdj.layout.consensusSequenceString().substring(querySeqRange))
            blastnRunDataMap[blastnRunData.key] = blastnRunData
        }

        // run blastn on those
        val blastnResults = BlastnUtil.runBlastn(
            sampleId, blastDir, blastDb,
            blastnRunDataMap.mapValues { runData -> runData.value.querySeq },
            outputDir, numThreads, BLASTN_MAX_EVALUE)

        // put all into an identity hash multimap
        val vdjToBlastnMatch: Multimap<BlastnRunData, BlastnMatch> = Multimaps.newListMultimap(IdentityHashMap()) { ArrayList() }

        for ((vdjKey, match) in blastnResults.entries())
        {
            val blastnRunData = blastnRunDataMap[vdjKey]

            if (blastnRunData == null)
            {
                sLogger.fatal("error processing blastn results: cannot find key: {}", vdjKey)
                throw RuntimeException("error processing blastn results: cannot find key: $vdjKey")
            }

            vdjToBlastnMatch.put(blastnRunData, match)
        }

        val blastnAnnotations = processBlastnMatches(blastnRunDataMap.values, vdjToBlastnMatch)

        // write blastn match tsv
        BlastnMatchTsvWriter.write(outputDir, sampleId, blastnAnnotations)

        return blastnAnnotations
    }

    // process the blastn matches for each VDJ, and set the blastnAnnotation in the VdjAnnotation
    // NOTE: we cannot use blastnMatches.keySet(), as it might not include some VDJs that returned no match
    fun processBlastnMatches(blastnRunDataList: Collection<BlastnRunData>, blastnMatches: Multimap<BlastnRunData, BlastnMatch>)
    : Collection<BlastnAnnotation>
    {
        val blastnAnnotations = ArrayList<BlastnAnnotation>()
        for (blastnRunData in blastnRunDataList)
        {
            blastnAnnotations.add(processBlastnMatches(blastnRunData, blastnMatches[blastnRunData]))
        }
        return blastnAnnotations
    }

    //
    fun processBlastnMatches(blastnRunData: BlastnRunData, blastnMatches: Collection<BlastnMatch>)
    : BlastnAnnotation
    {
        val vdjSequence: VDJSequence = blastnRunData.vdj

        val alignStartOffset = blastnRunData.querySeqRange.start

        // the matches are grouped by contigs, and each contig the matches are sorted by expected value
        // however we want to process matches by expected value first regardless of contig, so we need to
        // resort the matches

        // we also need to fix up the matches, since we did not use the full sequence to query blastn, the queryAlignStart
        // and queryAlignEnd indices are off
        blastnMatches.forEach { blastnMatch ->
            blastnMatch.queryAlignStart = blastnMatch.queryAlignStart + alignStartOffset
            blastnMatch.queryAlignEnd = blastnMatch.queryAlignEnd + alignStartOffset
        }

        val sortedMatches: List<BlastnMatch> = blastnMatches.sortedBy { m -> m.expectedValue }

        // we freeze the locus here. Reason is that there are cases where a low identity match (92%) from another
        // locus supercedes a 100% identity match from the correct locus
        val locus: IgTcrLocus? = vdjSequence.vAnchor?.geneType?.locus ?: vdjSequence.jAnchor?.geneType?.locus

        require(locus != null)

        var vGene: IgTcrGene? = null
        var vMatch: BlastnMatch? = null
        var dGene: IgTcrGene? = null
        var dMatch: BlastnMatch? = null
        var jGene: IgTcrGene? = null
        var jMatch: BlastnMatch? = null

        for (blastnMatch in sortedMatches)
        {
            // check if it matches whole way. Require 95% sequence identity
            if (blastnMatch.percentageIdent >= CiderConstants.BLASTN_MATCH_FULL_MATCH_IDENTITY &&
                blastnMatch.querySeqLen <= (blastnMatch.queryAlignEnd - blastnMatch.queryAlignStart) + 5)
            {
                // sLogger.debug("blastn matches ref genome: {}", blastnMatch.subjectTitle)
                // sLogger.debug("  query seq: {}", blastnMatch.alignedPartOfQuerySeq)
                // sLogger.debug("subject seq: {}", blastnMatch.alignedPartOfSubjectSeq)

                return BlastnAnnotation(
                    vdjSequence = vdjSequence,
                    fullMatch = blastnMatch, blastnStatus = BlastnStatus.NO_REARRANGEMENT)
            }

            var vdjGene: IgTcrGene? = findGene(blastnMatch)

            // for V/J gene segments, we mandate 90% identity
            if (vdjGene != null && vdjGene.region in arrayOf(IgTcrRegion.V_REGION, IgTcrRegion.J_REGION) && blastnMatch.percentageIdent < CiderConstants.BLASTN_MATCH_MIN_VJ_IDENTITY)
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
                    if (vdjGene.region == IgTcrRegion.V_REGION && isBetterMatch(vMatch, vGene, blastnMatch, vdjGene))
                    {
                        vGene = vdjGene
                        vMatch = blastnMatch
                    }
                    if (vdjGene.region == IgTcrRegion.D_REGION && isBetterMatch(dMatch, dGene, blastnMatch, vdjGene))
                    {
                        dGene = vdjGene
                        dMatch = blastnMatch
                    }
                    if (vdjGene.region == IgTcrRegion.J_REGION && isBetterMatch(jMatch, jGene, blastnMatch, vdjGene))
                    {
                        jGene = vdjGene
                        jMatch = blastnMatch
                    }
                }
            }
        }

        // determine status
        val blastnStatus: BlastnStatus = if (vGene != null)
        {
            if (dGene != null)
            {
                if (jGene != null)
                {
                    BlastnStatus.V_D_J
                }
                else
                {
                    BlastnStatus.V_D
                }
            }
            else if (jGene != null)
            {
                BlastnStatus.V_J
            }
            else
            {
                BlastnStatus.V_ONLY
            }
        }
        else
        {
            if (dGene != null)
            {
                if (jGene != null)
                {
                    BlastnStatus.D_J
                }
                else
                {
                    BlastnStatus.D_ONLY
                }
            }
            else if (jGene != null)
            {
                BlastnStatus.J_ONLY
            }
            else
            {
                BlastnStatus.NO_VDJ_ALIGNMENT
            }
        }

        return BlastnAnnotation(
            vdjSequence = vdjSequence,
            vGene = vGene,
            dGene = dGene,
            jGene = jGene,
            vMatch = vMatch,
            dMatch = dMatch,
            jMatch = jMatch,
            blastnStatus = blastnStatus)
    }

    fun findGene(blastnMatch: BlastnMatch) : IgTcrGene?
    {
        val matchLocation = toGenomicLocation(blastnMatch) ?: return null

        val chromosome = matchLocation.chromosome

        val geneDataList = vdjGenes[Pair(chromosome, blastnMatch.subjectFrame)]

        var bestGene : IgTcrGene? = null

        for (gene in geneDataList)
        {
            val geneLocation = gene.geneLocation ?: continue

            require(geneLocation.chromosome == chromosome)
            require(geneLocation.strand == blastnMatch.subjectFrame)

            if (bestGene == null || !bestGene.isFunctional)
            {
                // check if they overlap. We prioritise functional genes
                if (geneLocation.posStart <= matchLocation.posEnd + GENE_REGION_TOLERANCE &&
                    geneLocation.posEnd >= matchLocation.posStart - GENE_REGION_TOLERANCE)
                {
                    bestGene = gene
                }
            }
        }

        return bestGene
    }

    companion object
    {
        // from my test, it evalue of 1 can only match minimum 20 bases. If we want to match D segment that is shorter
        // we will need a higher cut off, maybe 10, but will get many false positive hits that are longer but more mismatches
        const val BLASTN_MAX_EVALUE = 1.0

        const val FLANKING_BASES = 50

        // D genes often are very short, for example, TRBD1 is only 12 bases. We allow more leeway to match
        // an alignment to the gene
        const val GENE_REGION_TOLERANCE = 50

        private val sLogger = LogManager.getLogger(BlastnAnnotator::class.java)

        fun blastnQuerySeqRange(vdj: VDJSequence) : IntRange
        {
            val fullSeq = vdj.layout.consensusSequenceString()
            val vdjSeq = vdj.sequence
            return blastnQuerySeqRange(vdjSeq, fullSeq)
        }

        fun blastnQuerySeqRange(vdjSeq: String, fullSeq: String): IntRange
        {
            val i = fullSeq.indexOf(vdjSeq)
            require(i >= 0)
            val range = Math.max(i - FLANKING_BASES, 0) until
                    Math.min(i + vdjSeq.length + FLANKING_BASES, fullSeq.length)
            return range
        }

        fun isBetterMatch(existingMatch: BlastnMatch?, existingGene: IgTcrGene?, newMatch: BlastnMatch, newGene: IgTcrGene) : Boolean
        {
            if (existingMatch == null || existingGene == null)
            {
                return true
            }

            if (Doubles.greaterThan(newMatch.bitScore, existingMatch.bitScore))
            {
                // always prefer higher bit score
                return true
            }

            if (Doubles.equal(newMatch.bitScore, existingMatch.bitScore) &&
                !existingGene.isFunctional &&
                newGene.isFunctional)
            {
                // if bitscores are equal, we prefer the functional one
                sLogger.trace("prefer functional gene: {}, bitscore: {} over non functional: {}, bitscore: {}",
                    newGene.geneAllele, newMatch.bitScore, existingGene.geneAllele, existingMatch.bitScore)
                return true
            }

            return false
        }
    }
}
