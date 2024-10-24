package com.hartwig.hmftools.cider.blastn

import com.google.common.collect.Multimap
import com.hartwig.hmftools.cider.CiderConstants
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.blastn.BlastnRunner
import com.hartwig.hmftools.common.genome.region.Strand

object BlastnUtil
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
            .build()
            .run(sequences)
    }
}
