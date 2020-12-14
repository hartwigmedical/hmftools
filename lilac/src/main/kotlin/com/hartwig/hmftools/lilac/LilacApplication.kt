package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.lilac.kmer.BamKmer
import com.hartwig.hmftools.lilac.kmer.HlaKmer
import com.hartwig.hmftools.lilac.prot.ProteinSequenceBoundaries
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile.firstFourDigits
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile.inflate
import htsjdk.samtools.QueryInterval
import htsjdk.samtools.SamReaderFactory
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.concurrent.Executors
import java.util.concurrent.Future

fun main(args: Array<String>) {
    LilacApplication().use { x -> x.run() }
}


class LilacApplication : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis()

    val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("LILAC-%d").build()
    val executorService = Executors.newFixedThreadPool(7, namedThreadFactory)


    override fun run() {

        val kmerLength = 10
        val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
        val sequenceAFile = "${resourcesDir}/A_prot.txt"
        val sequenceBFile = "${resourcesDir}/B_prot.txt"
        val sequenceCFile = "${resourcesDir}/C_prot.txt"
        val bamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.bam"
        // 6:29812305-33068473

        writeUnwrappedFiles(sequenceAFile, sequenceBFile, sequenceCFile)

        logger.info("Reading protein alignment files")
        val aSequences = ProteinSequenceFile.readFile(sequenceAFile).inflate().firstFourDigits()
        val bSequences = ProteinSequenceFile.readFile(sequenceBFile).inflate().firstFourDigits()
        val cSequences = ProteinSequenceFile.readFile(sequenceCFile).inflate().firstFourDigits()


        logger.info("Extracting ${kmerLength} length kmers")
        val aKmerMap = aSequences.map { Pair(it, it.uniqueExonicKmers(kmerLength, ProteinSequenceBoundaries.aBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()
        val bKmerMap = bSequences.map { Pair(it, it.uniqueExonicKmers(kmerLength, ProteinSequenceBoundaries.bBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()
        val cKmerMap = cSequences.map { Pair(it, it.uniqueExonicKmers(kmerLength, ProteinSequenceBoundaries.cBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()

        val hlaKmers = HlaKmer(aKmerMap, bKmerMap, cKmerMap)
        println(hlaKmers.uniqueKmers().size)
        println(hlaKmers.kmers().size)

//        val specificSequences = setOf(HlaAllele("A*01:01:01:01"), HlaAllele("A*11:01:01:01"))
//        val specificSequences = setOf(HlaAllele("B*08:01:01:01"), HlaAllele("B*56:01:01:01"))
//        val specificSequences = setOf(HlaAllele("C*01:02:01:01"), HlaAllele("C*07:01:01:01"))

        val kmerCount = kmerCounts(bamFile, hlaKmers.kmers())
        println(kmerCount.kmerCount().size)
        println(kmerCount.kmerCount().values.toIntArray().sum())


        val bamKmers = kmerCount.kmerCount().filter { it.value > 1 }.keys
        for ((proteinSequence, kmers) in hlaKmers.sequences()) {

//            logger.info("Processing ${proteinSequence.allele.fourDigitName()}")
//            for (kmer in kmers) {
//                logger.info("  Kmer $kmer count -> ${kmerCount.kmerCount()[kmer]}")
//            }


            if (kmers.isNotEmpty() && bamKmers.containsAll(kmers)) {
                val matches = kmerCount.kmerCount().filter { it.key in kmers }.map { it.value }.sum()
//                if (matches > 15000) {
                    logger.info("Found potential candidate ${proteinSequence.allele.toString()} -> $matches")
//                }

//                logger.info("PROCESSING ${proteinSequence.allele.fourDigitName()}")
//                for (kmer in kmers) {
//                    logger.info("       Kmer $kmer count -> ${kmerCount.kmerCount()[kmer]}")
//                }
//
            }


        }

        val uniqueKmers = bamKmers intersect hlaKmers.uniqueKmers()
        println(uniqueKmers.size)

        for (kmer in uniqueKmers) {
            val matches = kmerCount.kmerCount().filter { it.key in kmer }.map { it.value }.sum()
            logger.info("Found unique candidate ${hlaKmers.proteinSequence(kmer).allele.toString()} -> $matches -> $kmer")
        }


    }

    private fun kmerCounts(bamFile: String, kmers: Set<String>): BamKmer {
        val kmerCount = BamKmer(kmers)
        val futures = mutableListOf<Future<*>>()

        SamReaderFactory.makeDefault().open(File(bamFile)).use {

            val sequenceIndex = it.fileHeader.sequenceDictionary.getSequence("6").sequenceIndex
            val hlaa = QueryInterval(sequenceIndex, 29909331, 29914232)
            val hlac = QueryInterval(sequenceIndex, 31235946, 31240848)
            val hlab = QueryInterval(sequenceIndex, 31321260, 31325935)

//            val everything = QueryInterval(sequenceIndex, 29910331, 31324935)
//            val hlaa = QueryInterval(sequenceIndex, 29910331, 29913232)
//            val hlac = QueryInterval(sequenceIndex, 31236946, 31239848)
//            val hlab = QueryInterval(sequenceIndex, 31322260, 31324935)

            logger.info("Reading BAM file $bamFile")
            for (samRecord in it.queryOverlapping(arrayOf(hlaa, hlac, hlab))) {
                futures.add(executorService.submit { kmerCount.accept(samRecord) })
            }
        }

        logger.info("Searching for kmers")
        futures.forEach { it.get() }
        return kmerCount
    }


    private fun writeUnwrappedFiles(aFile: String, bFile: String, cFile: String) {
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.aBoundaries(), aFile, aFile.replace("txt", "unwrapped.txt"))
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.bBoundaries(), bFile, bFile.replace("txt", "unwrapped.txt"))
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.cBoundaries(), cFile, cFile.replace("txt", "unwrapped.txt"))
    }


    override fun close() {
        executorService.shutdown()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

}