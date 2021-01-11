package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.utils.SuffixTree
import com.hartwig.hmftools.lilac.coverage.ProteinCoverage
import com.hartwig.hmftools.lilac.ext.rollingKmers
import com.hartwig.hmftools.lilac.ext.toAminoAcids
import com.hartwig.hmftools.lilac.prot.ProteinSequence
import com.hartwig.hmftools.lilac.prot.ProteinSequenceBoundaries
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile.firstFourDigits
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile.inflate
import htsjdk.samtools.QueryInterval
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.Future

fun main(args: Array<String>) {
    LilacBamFilter().use { x -> x.run() }
}


class LilacBamFilter : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis()

    val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("LILAC-%d").build()
    val executorService = Executors.newFixedThreadPool(7, namedThreadFactory)


    fun coverage(boundaries: List<Int>, proteins: List<ProteinSequence>): List<ProteinCoverage> {
        return proteins.map { ProteinCoverage(it.allele, it.exonicProteins(boundaries).filter { it.length >= 40 }) }
    }

    override fun run() {

        // /data/common/tools/sambamba_v0.6.5/sambamba slice -o /data/experiments/201210_hla/GIABvsSELFv004R.hla.bam /data/data_archive/validation_data/pipeline_runs/200910_GIABvsSELFv004_v5.13/GIABvsSELFv004R/aligner/GIABvsSELFv004R.bam 6:29812305-33068473

        val minKmerLength = 10
        val maxKmerLength = 10
        val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
        val sequenceAFile = "${resourcesDir}/A_prot.txt"
        val sequenceBFile = "${resourcesDir}/B_prot.txt"
        val sequenceCFile = "${resourcesDir}/C_prot.txt"
        val bamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.bam"
        val outputBamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.filtered2.bam"

        logger.info("Reading protein alignment files")
        val aSequences = ProteinSequenceFile.readFile(sequenceAFile).inflate().firstFourDigits()
        val bSequences = ProteinSequenceFile.readFile(sequenceBFile).inflate().firstFourDigits()
        val cSequences = ProteinSequenceFile.readFile(sequenceCFile).inflate().firstFourDigits()



        logger.info("Extracting kmers between ${minKmerLength} and ${maxKmerLength} in length")
        val aKmerMap = aSequences.map { Pair(it, it.exonicKmers(minKmerLength, maxKmerLength, ProteinSequenceBoundaries.aBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()
        val bKmerMap = bSequences.map { Pair(it, it.exonicKmers(minKmerLength, maxKmerLength, ProteinSequenceBoundaries.bBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()
        val cKmerMap = cSequences.map { Pair(it, it.exonicKmers(minKmerLength, maxKmerLength, ProteinSequenceBoundaries.cBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()

        val allKmersSet = mutableSetOf<String>()
        aKmerMap.values.forEach { allKmersSet.addAll(it) }
        bKmerMap.values.forEach { allKmersSet.addAll(it) }
        cKmerMap.values.forEach { allKmersSet.addAll(it) }


        logger.info("Constructing suffix trees")
        val trees = allKmersSet.map { SuffixTree(it) }
        println(trees.size)

        val recordsToKeep = recordsToKeep(bamFile, trees)
        writeBam(bamFile, outputBamFile, recordsToKeep)

        println(recordsToKeep.size)

    }

    private fun writeBam(readBamFile: String, writeBamFile: String, keepers: Set<String>) {

        val input = SamReaderFactory.makeDefault().open(File(readBamFile))
        val output = SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(input.fileHeader, true, File(writeBamFile))

        val sequenceIndex = input.fileHeader.sequenceDictionary.getSequence("6").sequenceIndex
        val hlaa = QueryInterval(sequenceIndex, 29812305, 33068473)

        for (samRecord in input.queryOverlapping(arrayOf(hlaa))) {
            if (keepers.contains(samRecord.readName) || keepers.contains(samRecord.pairedReadName)) {
                output.addAlignment(samRecord)
            }
        }

        input.close()
        output.close()
    }

    private fun recordsToKeep(bamFile: String, trees: Collection<SuffixTree>): Set<String> {
        val futures = mutableListOf<Future<String?>>()
        SamReaderFactory.makeDefault().open(File(bamFile)).use {
            val sequenceIndex = it.fileHeader.sequenceDictionary.getSequence("6").sequenceIndex
            val hlaa = QueryInterval(sequenceIndex, 29812305, 33068473)

            LilacApplication.logger.info("Reading BAM file $bamFile")
            for (samRecord in it.queryOverlapping(arrayOf(hlaa))) {
                futures.add(executorService.submit(callableKeep(samRecord, trees)))
            }
        }

        LilacApplication.logger.info("Searching for kmers")
        return futures.map { it.get() }.filterNotNull().toSet()
    }


    fun callableKeep(record: SAMRecord, trees: Collection<SuffixTree>): Callable<String?> {
        return Callable {
            if (keep(record, trees)) {
//                logger.info("Processing YES ${record.readName}")
                return@Callable record.readName
            } else {
//                logger.info("Processing NO ${record.readName}")
                return@Callable null
            }
        }
    }


    fun keep(record: SAMRecord, trees: Collection<SuffixTree>): Boolean {
        val aminoAcids = record.toAminoAcids()
        val kmers = aminoAcids.flatMap { it.rollingKmers(10) }.toSet()

        for (kmer in kmers) {
            if (trees.any { it.contains(kmer) }) {
                return true
            }
        }
        return false
    }


    override fun close() {
        executorService.shutdown()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

}