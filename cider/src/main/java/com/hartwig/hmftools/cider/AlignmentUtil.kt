package com.hartwig.hmftools.cider

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.blastn.BlastnRunner
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
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

    // This limits the number of hit each query can get. Necessary to protect against edge cases
    const val MAX_TARGET_SEQUENCES = 5000

    // val REF_GENOME_NAME = "GRCh38.p13"
    val BLASTN_CHROMOSOME_ASSEMBLY_REGEX = Regex("^Homo sapiens chromosome (\\w+).*, GRCh38.p13 (.+)$")
    val BLASTN_PRIMARY_ASSEMBLY_NAME = "Primary Assembly".intern()

    fun toGenomicLocation(blastnMatch: BlastnMatch) : GenomicLocation?
    {
        // parse the chromosome / assembly
        val m = BLASTN_CHROMOSOME_ASSEMBLY_REGEX.matchEntire(blastnMatch.subjectTitle)

        if (m == null)
        {
            // might not have a chromosome, or it is mitochondrion
            return null
        }

        val chromosome = CiderConstants.BLAST_REF_GENOME_VERSION.versionedChromosome(m.groupValues[1]).intern()
        val assembly = m.groupValues[2].intern()

        val start: Int
        val end: Int

        when (blastnMatch.subjectFrame)
        {
            Strand.FORWARD -> { start = blastnMatch.subjectAlignStart; end = blastnMatch.subjectAlignEnd }
            Strand.REVERSE -> { start = blastnMatch.subjectAlignEnd; end = blastnMatch.subjectAlignStart }
        }

        val altAssemblyName = if (assembly == BLASTN_PRIMARY_ASSEMBLY_NAME) null else assembly

        return GenomicLocation(chromosome, start, end, blastnMatch.subjectFrame, altAssemblyName)
    }

    fun runBlastn(sampleId: String, blastDir: String, blastDb: String, sequences: Map<Int, String>, outputDir: String, numThreads: Int,
                  expectedValueCutoff: Double)
            : Multimap<Int, BlastnMatch>
    {
        return BlastnRunner.Builder()
            .withTask("blastn")
            .withPrefix(sampleId)
            .withBlastDir(blastDir)
            .withBlastDb(blastDb)
            .withOutputDir(outputDir)
            .withWordSize(WORD_SIZE)
            .withReward(MATCH_SCORE)
            .withPenalty(MISMATCH_SCORE)
            .withGapOpen(-GAP_OPENING_SCORE)
            .withGapExtend(-GAP_EXTEND_SCORE)
            .withExpectedValueCutoff(expectedValueCutoff)
            .withNumThreads(numThreads)
            .withSubjectBestHit(true)
            .withMaxTargetSeqs(MAX_TARGET_SEQUENCES)
            .build()
            .run(sequences)
    }

    data class BwaMemAlignment(
        val querySeq: String,
        val queryAlignStart: Int,
        val queryAlignEnd: Int,
        val refContig: String,
        val refStart: Int,
        val refEnd: Int,
        val strand: Strand,
        val alignmentScore: Int,
        val percentageIdent: Double,
    ) {
        init {
            require(refStart <= refEnd)
            require(queryAlignStart <= queryAlignEnd)
        }
    }

    fun parseChromosome(contig: String): String
    {
        // Parses ref genome contig to the primary assembly chromosome. E.g.:
        // chr14_KI270726v1_random
        return contig.split("_")[0]
    }

    fun runBwaMem(sequences: Map<Int, String>, refGenomeDictPath: String, refGenomeIndexPath: String, alignScoreThreshold: Int, numThreads: Int):
            Multimap<Int, BwaMemAlignment>
    {
        val refGenSeqDict = ReferenceSequenceFileFactory.loadDictionary(FileInputStream(refGenomeDictPath))
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
        aligner.dropRatioOption = 0.25f
        //aligner.maxMemIntvOption = 1000
        //aligner.maxSeedOccurencesOption = 1000
        //aligner.splitFactorOption = 0.5f
        //aligner.maxChainGapOption = 1000
        //aligner.zDropOption = 1000

        val keys = sequences.keys.toList()
        val seqs = keys.map { k -> sequences[k]!!.toByteArray() }
        val alignments = aligner.alignSeqs(seqs)
        val results = ArrayListMultimap.create<Int, BwaMemAlignment>()
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
                val querySeq = seqs[key.index]
                // Apparently these need to be inverted if the alignment reports the reverse strand
                // Do this to match Blastn behaviour
                val queryAlignStart = if (strand == Strand.FORWARD) { alignment.seqStart } else { querySeq.size - alignment.seqEnd } + 1
                val queryAlignEnd = if (strand == Strand.FORWARD) { alignment.seqEnd } else { querySeq.size - alignment.seqStart }
                require(queryAlignStart <= queryAlignEnd)
                // nMismatches is not the best name - it's actually the edit distance.
                // Which means this calculation is correct for mismatches and gaps.
                val percentIdentity = 100 * (1 - (alignment.nMismatches.toDouble() / (queryAlignEnd - queryAlignStart + 1)))
                val resAlignment = BwaMemAlignment(
                    sequences[key.value]!!,
                    queryAlignStart,
                    queryAlignEnd,
                    refContig,
                    refStart,
                    refEnd,
                    strand,
                    alignment.alignerScore,
                    percentIdentity
                )
                results.put(key.value, resAlignment)
            }
        }
        return results
    }
}