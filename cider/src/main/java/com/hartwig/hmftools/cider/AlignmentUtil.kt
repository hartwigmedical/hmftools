package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.CiderConstants.ALIGNMENT_BATCH_SIZE
import com.hartwig.hmftools.cider.genes.Contig
import com.hartwig.hmftools.cider.genes.ContigType
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.region.BaseRegion
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import org.apache.logging.log4j.LogManager
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex
import java.io.FileInputStream

object AlignmentUtil
{
    // For the scoring function, the match/mismatch score 1/-4 optimizes the scoring for 100% identical sequences and 1/-1 for 75% identical
    // sequences. The default for NCBI Blastn is 2/-3, which is optimal for 89% identical sequences. BWA uses 1/-4.
    // There is also gap opening and gap extension. BWA uses gap opening of -6 and gap extension of -1.
    // For blastn default scoring, see: https://www.ncbi.nlm.nih.gov/books/NBK279684/

    // we set it to 1/-4/-5/-2 which is optimal for 100% identical sequences
    // My test shows that this scoring would mostly prefer shorter matches with higher identity than longer matches with lower identity.
    const val MATCH_SCORE = 1
    const val MISMATCH_SCORE = -4
    const val GAP_OPENING_SCORE = -5
    const val GAP_EXTEND_SCORE = -2

    const val WORD_SIZE = 9

    private val sLogger = LogManager.getLogger(AlignmentUtil::class.java)

    data class BwaMemAlignment(
        val querySeq: String,
        val queryAlignStart: Int,
        val queryAlignEnd: Int,
        val refContig: String,      // For GRCh38 this is prefixed with "chr", for GRCh37 it is not.
        val refStart: Int,
        val refEnd: Int,
        val refStrand: Strand,
        val alignmentScore: Int,
        val editDistance: Int,
        val percentageIdent: Double,
    ) {
        init {
            require(refStart <= refEnd)
            require(queryAlignStart <= queryAlignEnd)
        }
    }

    fun toGenomicLocation(alignment: BwaMemAlignment): GenomicLocation?
    {
        val contig = Contig(alignment.refContig)

        if (contig.type == ContigType.MITOCHONDRIAL || contig.type == ContigType.UNPLACED)
        {
            // This matches previous behaviour with Blast.
            // TODO: review?
            return null
        }

        return GenomicLocation(contig, BaseRegion(alignment.refStart, alignment.refEnd), alignment.refStrand)
    }

    fun runBwaMem(sequences: Map<Int, String>, refGenomeDictPath: String, refGenomeIndexPath: String, alignScoreThreshold: Int, numThreads: Int):
            Map<Int, ArrayList<BwaMemAlignment>>
    {
        sLogger.debug("Aligning ${sequences.size} sequences")

        val refGenSeqDict = ReferenceSequenceFileFactory.loadDictionary(FileInputStream(refGenomeDictPath))
        val aligner = createBwaMemAligner(refGenomeIndexPath, alignScoreThreshold, numThreads)

        val keys = sequences.keys.sorted()
        val results = HashMap<Int, ArrayList<BwaMemAlignment>>()
        // Alignments are batches because with our BWA-MEM settings, too much memory is allocated with large BAMs.
        for (i in 0 until keys.size step ALIGNMENT_BATCH_SIZE) {
            val batchKeys = keys.subList(i, minOf(i + ALIGNMENT_BATCH_SIZE, keys.size))
            sLogger.debug("Running BWA-MEM on batch of ${batchKeys.size} sequences")
            val batchByteSeqs = batchKeys.map { k -> sequences[k]!!.toByteArray() }
            // Perform a JVM garbage collection before alignment to reduce memory pressure.
            System.gc()
            val batchAlignments = aligner.alignSeqs(batchByteSeqs)
            val batchResults = parseBwaMemAlignments(sequences, batchKeys, batchAlignments, refGenSeqDict)
            batchResults.forEach { (key, value) -> results.computeIfAbsent(key) { ArrayList() }.addAll(value) }
        }
        sLogger.debug("Alignment complete")
        return results
    }

    private fun createBwaMemAligner(refGenomeIndexPath: String, alignScoreThreshold: Int, numThreads: Int): BwaMemAligner
    {
        val index = BwaMemIndex(refGenomeIndexPath)
        val aligner = BwaMemAligner(index)
        aligner.nThreadsOption = numThreads
        aligner.flagOption = aligner.flagOption or BwaMemAligner.MEM_F_ALL
        aligner.outputScoreThresholdOption = alignScoreThreshold
        aligner.minSeedLengthOption = WORD_SIZE
        aligner.matchScoreOption = MATCH_SCORE
        aligner.mismatchPenaltyOption = -MISMATCH_SCORE
        aligner.dGapOpenPenaltyOption = -GAP_OPENING_SCORE
        aligner.iGapOpenPenaltyOption = -GAP_OPENING_SCORE
        aligner.dGapExtendPenaltyOption = -GAP_EXTEND_SCORE
        aligner.iGapExtendPenaltyOption = -GAP_EXTEND_SCORE
        // Relax pruning parameters to encourage more alignments to be found.
        // Otherwise in some cases BWA will miss alignments which we know are correct for gene annotation.
        // Note these parameters reduce performance and cause large transient memory allocations.
        aligner.dropRatioOption = 0.1f
        aligner.splitFactorOption = 0.5f
        aligner.maxMemIntvOption = 500
        aligner.maxSeedOccurencesOption = 2000
        aligner.maxChainGapOption = 300
        aligner.zDropOption = 300
        return aligner
    }

    private fun parseBwaMemAlignments(sequences: Map<Int, String>, keys: List<Int>,
                                      alignments: List<List<org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment>>,
                                      refGenSeqDict: SAMSequenceDictionary): Map<Int, List<BwaMemAlignment>>
    {
        val results = HashMap<Int, ArrayList<BwaMemAlignment>>()
        for (key in keys.withIndex()) {
            for (alignment in alignments[key.index]) {
                if (alignment.samFlag and 0x4 != 0)
                {
                    // Not a real alignment, means the query is unmapped.
                    continue
                }

                val refContig = refGenSeqDict.getSequence(alignment.refId).sequenceName
                val refStart = alignment.refStart + 1   // apparently BWA lib gives 0-based index
                val refEnd = alignment.refEnd
                require(refStart <= refEnd)
                val strand = if ((alignment.samFlag and 0x10) == 0) { Strand.FORWARD } else { Strand.REVERSE }
                val querySeq = sequences[key.value]!!
                // Apparently these need to be inverted if the alignment reports the reverse strand
                // Do this to match Blastn behaviour
                val queryAlignStart = if (strand == Strand.FORWARD) { alignment.seqStart } else { querySeq.length - alignment.seqEnd } + 1
                val queryAlignEnd = if (strand == Strand.FORWARD) { alignment.seqEnd } else { querySeq.length - alignment.seqStart }
                require(queryAlignStart <= queryAlignEnd)
                // nMismatches is not the best name - it's actually the edit distance.
                // Which means this calculation is correct for mismatches and gaps.
                val percentIdentity = 100 * (1 - (alignment.nMismatches.toDouble() / (queryAlignEnd - queryAlignStart + 1)))
                val resAlignment = BwaMemAlignment(
                    querySeq,
                    queryAlignStart,
                    queryAlignEnd,
                    refContig,
                    refStart,
                    refEnd,
                    strand,
                    alignment.alignerScore,
                    alignment.nMismatches,
                    percentIdentity
                )
                results.computeIfAbsent(key.value) { ArrayList() }.add(resAlignment)
            }
        }
        return results
    }
}