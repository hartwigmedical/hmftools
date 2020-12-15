package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.lilac.kmer.KmerCount
import com.hartwig.hmftools.lilac.prot.ProteinSequenceBoundaries
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile
import htsjdk.samtools.QueryInterval
import htsjdk.samtools.SamReaderFactory
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit

fun main(args: Array<String>) {
    LilacApplication().use { x -> x.run() }
}


class LilacApplication : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis()

    val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("SAGE-%d").build()
    val executorService = Executors.newFixedThreadPool(7, namedThreadFactory)


    override fun run() {

        val kmerLength = 10
        val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
        val sequenceAFile = "${resourcesDir}/A_prot.txt"
        val sequenceBFile = "${resourcesDir}/B_prot.txt"
        val sequenceCFile = "${resourcesDir}/C_prot.txt"

        val bamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.bam"
        // 6:29812305-33068473


        logger.info("Reading protein alignment files")
        val aSequences = ProteinSequenceFile.readWrappedFile("${resourcesDir}/A_prot.txt");
        val bSequences = ProteinSequenceFile.readWrappedFile("${resourcesDir}/B_prot.txt");
        val cSequences = ProteinSequenceFile.readWrappedFile("${resourcesDir}/C_prot.txt");

        logger.info("Extracting ${kmerLength} length kmers")
        val aSequenceKmers = aSequences.flatMap { it.exonicKmers(kmerLength, ProteinSequenceBoundaries.aBoundaries()) }
        val bSequenceKmers = bSequences.flatMap { it.exonicKmers(kmerLength, ProteinSequenceBoundaries.bBoundaries()) }
        val cSequenceKmers = cSequences.flatMap { it.exonicKmers(kmerLength, ProteinSequenceBoundaries.cBoundaries()) }

        val allKmers = mutableSetOf<String>()
        allKmers.addAll(aSequenceKmers)
        allKmers.addAll(bSequenceKmers)
        allKmers.addAll(cSequenceKmers)

        val kmerCount = KmerCount(allKmers)


        var i = 0;
        SamReaderFactory.makeDefault().open(File(bamFile)).use {

            val sequenceIndex = it.fileHeader.sequenceDictionary.getSequence("6").sequenceIndex

//            val everything = QueryInterval(sequenceIndex, 29910331, 31324935)
//            val hlaa = QueryInterval(sequenceIndex, 29910331, 29913232)
//            val hlac = QueryInterval(sequenceIndex, 31236946, 31239848)
//            val hlab = QueryInterval(sequenceIndex, 31322260, 31324935)

            val hlaa = QueryInterval(sequenceIndex, 29909331, 29914232)
            val hlac = QueryInterval(sequenceIndex, 31235946, 31240848)
            val hlab = QueryInterval(sequenceIndex, 31321260, 31325935)

            logger.info("Reading BAM")
            for (samRecord in it.queryOverlapping(arrayOf(hlaa, hlac, hlab))) {
                executorService.submit { kmerCount.accept(samRecord) }
//                kmerCount.accept(samRecord)
                i++
            }
        }

        logger.info("Processing codon matches")
        executorService.shutdown()
        executorService.awaitTermination(20, TimeUnit.MINUTES)

        println(i)
        println(allKmers.size)
        println(kmerCount.kmerCount().size)
        println(kmerCount.kmerCount().values.toIntArray().sum())

    }


    private fun writeUnwrappedFiles(aFile: String, bFile: String, cFile: String) {
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.aBoundaries(), aFile, aFile.replace("txt", "unwrapped.txt"))
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.bBoundaries(), bFile, bFile.replace("txt", "unwrapped.txt"))
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.cBoundaries(), cFile, cFile.replace("txt", "unwrapped.txt"))
    }


    override fun close() {
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

}