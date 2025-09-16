package com.hartwig.hmftools.cider.alignment

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.cider.CiderConstants
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.blastn.BlastnRunner
import com.hartwig.hmftools.common.codon.Nucleotides
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome
import com.hartwig.hmftools.common.genome.region.Strand
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex
import kotlin.math.max
import kotlin.math.min

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
    val CHROMOSOME_ASSEMBLY_REGEX = Regex("^Homo sapiens chromosome (\\w+).*, GRCh38.p13 (.+)$")
    val PRIMARY_ASSEMBLY_NAME = "Primary Assembly".intern()

    fun isPrimaryAssembly(assemblyName: String) : Boolean
    {
        return assemblyName == PRIMARY_ASSEMBLY_NAME
    }

    fun toGenomicLocation(blastnMatch: BlastnMatch) : GenomicLocation?
    {
        // parse the chromosome / assembly
        val m = CHROMOSOME_ASSEMBLY_REGEX.matchEntire(blastnMatch.subjectTitle)

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

        return GenomicLocation(chromosome, start, end, blastnMatch.subjectFrame, if (isPrimaryAssembly(assembly)) null else assembly)
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

    data class BwaMemMatch(
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

    fun runBwaMem(sequences: Map<Int, String>, refGenomeFastaPath: String, refGenomeIndexPath: String, numThreads: Int):
            Multimap<Int, BwaMemMatch>
    {
        val refGenome = loadRefGenome(refGenomeFastaPath)
        val refGenomeVersion = RefGenomeSource.deriveRefGenomeVersion(refGenome)
        val index = BwaMemIndex(refGenomeIndexPath)
        val aligner = BwaMemAligner(index)
        aligner.setNThreadsOption(numThreads)
        aligner.setFlagOption(aligner.getFlagOption() or BwaMemAligner.MEM_F_ALL)
        aligner.setOutputScoreThresholdOption(1)
        aligner.setMinSeedLengthOption(WORD_SIZE)
        aligner.setMatchScoreOption(MATCH_SCORE)
        aligner.setMismatchPenaltyOption(-MISMATCH_SCORE)
        aligner.dGapOpenPenaltyOption = -GAP_OPENING_SCORE;
        aligner.iGapOpenPenaltyOption = -GAP_OPENING_SCORE;
        aligner.dGapExtendPenaltyOption = -GAP_EXTEND_SCORE;
        aligner.iGapExtendPenaltyOption = -GAP_EXTEND_SCORE;
//        aligner.setMaxMemIntvOption(2000)
//        aligner.setMaxSeedOccurencesOption(2000)

        val keys = sequences.keys.toList()
        val seqs = keys.map { k -> sequences[k]!!.toByteArray() }
        val alignments = aligner.alignSeqs(seqs)
        val results = ArrayListMultimap.create<Int, BwaMemMatch>()
        for (key in keys.withIndex()) {
            for (alignment in alignments[key.index]) {
                val chromosome = if (alignment.refId >= 0 && alignment.refId < HumanChromosome.entries.size) {
                    refGenomeVersion.versionedChromosome(HumanChromosome.entries[alignment.refId].toString())
                }
                else {
                    continue
                }
                val refStart = alignment.refStart + 1   // apparently BWA lib gives 0-based index
                val refEnd = alignment.refEnd
                require(refStart <= refEnd)
                val strand = if ((alignment.samFlag and 0x10) == 0) { Strand.FORWARD } else { Strand.REVERSE }
                var refSeq = refGenome.getBases(chromosome, refStart, refEnd)
                refSeq = if (strand == Strand.FORWARD) { refSeq } else {
                    Nucleotides.reverseComplementBases(refSeq)
                }
                val querySeq = seqs[key.index]
                // Apparently these need to be inverted if the alignment reports the reverse strand
                // Do this to match Blastn behaviour
                val queryAlignStart = if (strand == Strand.FORWARD) { alignment.seqStart } else { querySeq.size - alignment.seqEnd }
                val queryAlignEnd = if (strand == Strand.FORWARD) { alignment.seqEnd } else { querySeq.size - alignment.seqStart }
                require(queryAlignStart <= queryAlignEnd)
                val alignedQuerySeq = querySeq.copyOfRange(queryAlignStart, queryAlignEnd)
                val percentIdentity = calcPercentIdentity(refSeq, alignedQuerySeq)
                val resAlignment = BwaMemMatch(
                    sequences[key.value]!!,
                    // Convert to 1-indexed for the rest of the code
                    queryAlignStart + 1,
                    queryAlignEnd,
                    chromosome,
                    alignment.refStart,
                    alignment.refEnd,
                    strand,
                    alignment.alignerScore,
                    percentIdentity
                )
                results.put(key.value, resAlignment)
            }
        }
        return results
    }

    // TODO: this is wrong because it doesn't take into account gaps. Probably need to parse CIGAR score instead.
    private fun calcPercentIdentity(refSeq: ByteArray, alignedQuerySeq: ByteArray): Double {
        var match = 0
        for (i in 0 until min(refSeq.size, alignedQuerySeq.size)) {
            if (refSeq[i] == alignedQuerySeq[i]) {
                ++match
            }
        }
        return 100 * match.toDouble() / max(refSeq.size, alignedQuerySeq.size)
    }
}
