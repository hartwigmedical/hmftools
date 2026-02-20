package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.CiderConstants.ALIGNER_GAP_EXTEND_SCORE
import com.hartwig.hmftools.cider.CiderConstants.ALIGNER_GAP_OPENING_SCORE
import com.hartwig.hmftools.cider.CiderConstants.ALIGNER_MATCH_SCORE
import com.hartwig.hmftools.cider.CiderConstants.ALIGNER_MISMATCH_SCORE
import com.hartwig.hmftools.cider.CiderConstants.ALIGNER_WORD_SIZE
import com.hartwig.hmftools.cider.CiderConstants.BLASTN_CHROMOSOME_ASSEMBLY_REGEX
import com.hartwig.hmftools.cider.CiderConstants.BLASTN_MAX_TARGET_SEQUENCES
import com.hartwig.hmftools.cider.CiderConstants.BLASTN_PRIMARY_ASSEMBLY_NAME
import com.hartwig.hmftools.cider.CiderConstants.BWAMEM_BATCH_SIZE
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.blastn.BlastnRunner
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.CigarElement
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import org.apache.logging.log4j.LogManager
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex
import java.io.FileInputStream
import java.io.InputStream

private val sLogger = LogManager.getLogger("AlignmentUtil")

data class Alignment(
    val querySeq: String,   // Original query sequence input into the aligner.
    val queryRange: IntRange,  // 0-based, with respect to querySeq.
    val refContig: String,
    val refRange: IntRange, // 0-based.
    val strand: Strand,     // If REVERSE, then the reverse complement of querySeq was matched to the reference.
    val alignmentScore: Int,
    val editDistance: Int,  // With respect to querySeq.
    val cigar: List<CigarElement>,
    val refContigLength: Int,
)
{
    init
    {
        require(queryRange.start >= 0)
        require(queryRange.endInclusive < querySeq.length)
        require(queryRange.start <= queryRange.endInclusive)
        require(refRange.start >= 0)
        require(refRange.start <= refRange.endInclusive)
        require(refRange.endInclusive < refContigLength)
        require(editDistance >= startClip + endClip)
        require(cigar.isNotEmpty())
    }

    val queryAlignLength = queryRange.endInclusive + 1 - queryRange.start

    // Number of bases clipped at the start of querySeq.
    val startClip: Int get() = queryRange.start

    // Number of bases clipped at the end of querySeq.
    val endClip: Int get() = querySeq.length - (queryRange.endInclusive + 1)

    // Edit distance of the aligned subsequence of querySeq. I.e. excluding clipping.
    val alignedEditDistance: Int get() = editDistance - startClip - endClip
}

fun blastnMatchtoGenomicLocation(blastnMatch: BlastnMatch) : GenomicLocation?
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

fun runBlastn(sampleId: String, blastDir: String, blastDb: String, sequences: List<String>, outputDir: String, numThreads: Int,
              expectedValueCutoff: Double)
        : List<Collection<BlastnMatch>>
{
    var key = 0
    val sequencesByKey = sequences.associateBy { key++ }

    val resultsMap = BlastnRunner.Builder()
        .withTask("blastn")
        .withPrefix(sampleId)
        .withBlastDir(blastDir)
        .withBlastDb(blastDb)
        .withOutputDir(outputDir)
        .withWordSize(ALIGNER_WORD_SIZE)
        .withReward(ALIGNER_MATCH_SCORE)
        .withPenalty(ALIGNER_MISMATCH_SCORE)
        .withGapOpen(-ALIGNER_GAP_OPENING_SCORE)
        .withGapExtend(-ALIGNER_GAP_EXTEND_SCORE)
        .withExpectedValueCutoff(expectedValueCutoff)
        .withNumThreads(numThreads)
        .withSubjectBestHit(false)
        .withMaxTargetSeqs(BLASTN_MAX_TARGET_SEQUENCES)
        .build()
        .run(sequencesByKey)

    return sequences.indices.map { key -> resultsMap.get(key) }
}

fun runBwaMem(sequences: List<String>, refGenomeDictPath: String, refGenomeIndexPath: String, alignScoreThreshold: Int, numThreads: Int):
        List<List<Alignment>> =
    runBwaMem(sequences, FileInputStream(refGenomeDictPath), refGenomeIndexPath, alignScoreThreshold, numThreads)

fun runBwaMem(sequences: List<String>, refGenomeDict: InputStream, refGenomeIndexPath: String, alignScoreThreshold: Int, numThreads: Int):
        List<List<Alignment>> =
    runBwaMem(sequences, ReferenceSequenceFileFactory.loadDictionary(refGenomeDict), refGenomeIndexPath, alignScoreThreshold, numThreads)

fun runBwaMem(
    sequences: List<String>, refGenomeDict: SAMSequenceDictionary, refGenomeIndexPath: String, alignScoreThreshold: Int, numThreads: Int
): List<List<Alignment>>
{
    sLogger.debug("Aligning ${sequences.size} sequences with BWA-MEM")

    val aligner = createBwaMemAligner(refGenomeIndexPath, alignScoreThreshold, numThreads)

    val results = ArrayList<List<Alignment>>()
    // Alignments are batches because with our BWA-MEM settings, too much memory is allocated with large BAMs.
    for (i in 0 until sequences.size step BWAMEM_BATCH_SIZE) {
        val batchSequences = sequences.subList(i, minOf(i + BWAMEM_BATCH_SIZE, sequences.size))
        sLogger.debug("Running BWA-MEM on batch of ${batchSequences.size} sequences")
        val batchByteSeqs = batchSequences.map { it.toByteArray() }
        // Perform a JVM garbage collection before alignment to reduce memory pressure.
        System.gc()
        val batchAlignments = aligner.alignSeqs(batchByteSeqs)
        val batchResults = parseBwaMemAlignments(batchSequences, batchAlignments, refGenomeDict)
        results.addAll(batchResults)
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
    aligner.minSeedLengthOption = ALIGNER_WORD_SIZE
    aligner.matchScoreOption = ALIGNER_MATCH_SCORE
    aligner.mismatchPenaltyOption = -ALIGNER_MISMATCH_SCORE
    aligner.dGapOpenPenaltyOption = -ALIGNER_GAP_OPENING_SCORE
    aligner.iGapOpenPenaltyOption = -ALIGNER_GAP_OPENING_SCORE
    aligner.dGapExtendPenaltyOption = -ALIGNER_GAP_EXTEND_SCORE
    aligner.iGapExtendPenaltyOption = -ALIGNER_GAP_EXTEND_SCORE
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

private fun parseBwaMemAlignments(sequences: List<String>,
                                  alignments: List<List<org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment>>,
                                  refGenSeqDict: SAMSequenceDictionary): List<List<Alignment>>
{
    require(sequences.size == alignments.size)
    return sequences.zip(alignments)
        .map { (querySeq, queryAlignments) ->
            queryAlignments.mapNotNull { parseBwaMemAlignment(querySeq, it, refGenSeqDict) } }
}

private fun parseBwaMemAlignment(
    querySeq: String, alignment: org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment, refGenSeqDict: SAMSequenceDictionary)
    : Alignment?
{
    if (alignment.samFlag and 0x4 != 0)
    {
        // Not a real alignment, means the query is unmapped.
        return null
    }

    val refContigSequence = refGenSeqDict.getSequence(alignment.refId)
    val refContig = refContigSequence.sequenceName
    // Apparently BWA lib gives 0-based index
    val refRange = alignment.refStart until alignment.refEnd

    val strand = if ((alignment.samFlag and 0x10) == 0) Strand.FORWARD else Strand.REVERSE

    val cigar = cigarElementsFromStr(alignment.cigar)
    val (startClip, endClip) = getQueryClipping(cigar, strand)

    val queryRange = startClip until querySeq.length - endClip

    // nMismatches is not the best name - it's actually the edit distance but excluding clipping.
    val editDistance = alignment.nMismatches + startClip + endClip

    return Alignment(
        querySeq,
        queryRange,
        refContig,
        refRange,
        strand,
        alignment.alignerScore,
        editDistance,
        cigar,
        refContigSequence.sequenceLength
    )
}

private fun getQueryClipping(cigar: List<CigarElement>, strand: Strand): Pair<Int, Int>
{
    val leftClip = if (cigar[0].operator.isClipping) cigar[0].length else 0
    val rightClip = if (cigar.last().operator.isClipping) cigar.last().length else 0
    return if (strand == Strand.FORWARD) Pair(leftClip, rightClip) else Pair(rightClip, leftClip)
}
