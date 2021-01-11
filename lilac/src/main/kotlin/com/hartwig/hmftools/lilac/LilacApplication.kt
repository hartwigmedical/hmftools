package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.utils.SuffixTree
import com.hartwig.hmftools.lilac.coverage.BamCoverage
import com.hartwig.hmftools.lilac.coverage.BamCoverage3
import com.hartwig.hmftools.lilac.coverage.ProteinCoverage
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.kmer.BamKmer
import com.hartwig.hmftools.lilac.kmer.HlaKmer
import com.hartwig.hmftools.lilac.prot.ProteinSequence
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


    fun coverage(boundaries: List<Int>, proteins: List<ProteinSequence>): List<ProteinCoverage> {
        return proteins.map { ProteinCoverage(it.allele, it.exonicProteins(boundaries).filter { it.length >= 40 }) }
    }

    override fun run() {

        val minKmerLength = 10
        val maxKmerLength = 20
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

//        val jon = bSequences.filter { it.allele == HlaAllele("B*56:68") }[0]
//        println(jon.exonicProteins(ProteinSequenceBoundaries.bBoundaries()))

        for (i in 1 until aSequences.size) {
            println("$i of ${aSequences.size}")
            aSequences[i].exonicProteins(ProteinSequenceBoundaries.aBoundaries()).map { SuffixTree(it) }
        }



        logger.info("Extracting kmers between ${minKmerLength} and ${maxKmerLength} in length")
        val aKmerMap = aSequences.map { Pair(it, it.exonicKmers(minKmerLength, maxKmerLength, ProteinSequenceBoundaries.aBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()
        val bKmerMap = bSequences.map { Pair(it, it.exonicKmers(minKmerLength, maxKmerLength, ProteinSequenceBoundaries.bBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()
        val cKmerMap = cSequences.map { Pair(it, it.exonicKmers(minKmerLength, maxKmerLength, ProteinSequenceBoundaries.cBoundaries())) }.filter { it.second.isNotEmpty() }.toMap()

        val hlaKmers = HlaKmer(listOf(aKmerMap, bKmerMap, cKmerMap))


        for (i in hlaKmers.sequences) {
            for (j in hlaKmers.sequences) {
                if (i.key != j.key && i.value == j.value) {
                    logger.warn("${i.key.allele} is indistingusable from ${j.key.allele}")
                }
            }
        }


        logger.info("Kmers: ${hlaKmers.kmers().size}")
        logger.info("Unique kmers: ${hlaKmers.uniqueKmers().size}")

        val specificSequencesA = setOf(HlaAllele("A*01:01:01:01"), HlaAllele("A*11:01:01:01"))
        val specificSequencesB = setOf(HlaAllele("B*08:01:01:01"), HlaAllele("B*56:01:01:01"))
        val specificSequencesC = setOf(HlaAllele("C*01:02:01:01"), HlaAllele("C*07:01:01:01"))
        val actualSequences = mutableSetOf<HlaAllele>()
        actualSequences.addAll(specificSequencesA)
        actualSequences.addAll(specificSequencesB)
        actualSequences.addAll(specificSequencesC)

        val bComplexLong = setOf(
                HlaAllele("B*07:02:01:01"),
                HlaAllele("B*08:01:01:01"),
                HlaAllele("B*08:20:01"),
                HlaAllele("B*08:26:03"),
                HlaAllele("B*08:39"),
                HlaAllele("B*08:132"),
                HlaAllele("B*08:134"),
                HlaAllele("B*08:151"),
                HlaAllele("B*08:221"),
                HlaAllele("B*08:225"),
                HlaAllele("B*08:230"),
                HlaAllele("B*08:233"),
                HlaAllele("B*08:248"),
                HlaAllele("B*42:21"),
                HlaAllele("B*55:52"),
                HlaAllele("B*56:01:01:01"),
                HlaAllele("B*56:68"))


        val bComplex = setOf(HlaAllele("B*08:01:01:01"), HlaAllele("B*56:01:01:01"), HlaAllele("B*56:68"))
        val specificSequences = setOf(HlaAllele("C*01:02:01:01"), HlaAllele("C*07:01:01:01"), HlaAllele("C*07:877"), HlaAllele("C*07:879"), HlaAllele("C*07:882"))


        val kmerCount = kmerCounts(bamFile, hlaKmers.kmers())
        logger.info("Found ${kmerCount.size} distinct kmers with ${kmerCount.values.toIntArray().sum()} coverage")

        val initialCandidates = initialCandidates(hlaKmers, kmerCount)
        val initialCandidateAlleles = initialCandidates.map { it.key.allele }


        logger.info("*********** Coverage ***********")


        val aCoverage = coverage(ProteinSequenceBoundaries.aBoundaries(), aSequences).filter { it.hlaAllele in initialCandidateAlleles }
        val bCoverage = coverage(ProteinSequenceBoundaries.bBoundaries(), bSequences).filter { it.hlaAllele in initialCandidateAlleles }
        val cCoverage = coverage(ProteinSequenceBoundaries.cBoundaries(), cSequences).filter { it.hlaAllele in initialCandidateAlleles }
        val allCoverage = mutableListOf<ProteinCoverage>()
        allCoverage.addAll(aCoverage)
        allCoverage.addAll(bCoverage)
        allCoverage.addAll(cCoverage)



//        for (proteinCoverage in allCoverage) {
//            println(proteinCoverage.map.values.map { it.exonSequence })
//
//        }

        // Coverage
        calculateCoverage3(minKmerLength, bamFile, allCoverage)
        val secondaryCandidates = allCoverage.filter { it.isCovered() }

        println("${initialCandidates.size} initial candidates")
        println("${secondaryCandidates.size} secondary candidates")

        logger.info("Candidate\tMinCoverage\tMaxCoverage\tAvgCoverage")
        for (coverage in secondaryCandidates) {
            val minCoverage = coverage.minCoverage()
                logger.info("${coverage.hlaAllele.toString().padEnd(14, ' ')}\t${minCoverage}\t${coverage.maxCoverage()}\t${coverage.avgCoverage()}")

        }


        println("sdf")



//
//
//        logger.info("*********** Round 2 ***********")
//        val candidateKmers = HlaKmer(listOf(initialCandidates))
//        val uniqueCandidateKmers = candidateKmers.uniqueKmers
//        round2(candidateKmers, kmerCount.filter { it.key in uniqueCandidateKmers })
//
////        val uniqueKmers = bamKmers intersect hlaKmers.uniqueKmers()
////        println(uniqueKmers.size)
////
////        for (kmer in uniqueKmers) {
////            val matches = kmerCount.kmerCount().filter { it.key in kmer }.map { it.value }.sum()
////            logger.info("Found unique candidate ${hlaKmers.proteinSequence(kmer).allele.toString()} -> $matches -> $kmer")
////        }
//
//        logger.info("Found ${initialCandidates.size} candidates")

    }



    private fun initialCandidates(hlaKmers: HlaKmer, kmerCount: Map<String, Int>): Map<ProteinSequence, Set<String>> {
        val result = mutableMapOf<ProteinSequence, Set<String>>()

        val bamKmers = kmerCount.filter { it.value > 1 }.keys
        val uniqueBamKmers = bamKmers intersect hlaKmers.uniqueKmers


        logger.info("Candidate\tMinCoverage\tMaxCoverage\tTotalCoverage")
        for ((proteinSequence, kmers) in hlaKmers.sequences()) {
            val uniqueKmers = kmers intersect uniqueBamKmers

            if (kmers.isNotEmpty() && bamKmers.containsAll(kmers)) {
                result[proteinSequence] = kmers

                val totalMatches = kmerCount.filter { it.key in kmers }.map { it.value }.sum()
                val minMatch = kmerCount.filter { it.key in kmers }.map { it.value }.min()
                val maxMatch = kmerCount.filter { it.key in kmers }.map { it.value }.max()
                logger.info("${proteinSequence.allele.toString().padEnd(14, ' ')}\t$minMatch\t$maxMatch\t$totalMatches")
//                for (kmer in kmers) {
//                    logger.info("       Kmer $kmer count -> ${kmerCount[kmer]}")
//                }
//
                for (kmer in uniqueKmers) {
//                    logger.info("   includes unique kmer $kmer -> ${kmerCount[kmer]}")
                }
            }
        }

        return result
//                .filter { it.key.allele != HlaAllele("B*56:01:01:01") }
//                .filter { !it.key.allele.toString().endsWith("N") }
    }

    private fun round2(hlaKmers: HlaKmer, kmerCount: Map<String, Int>): Map<ProteinSequence, Set<String>> {
        val result = mutableMapOf<ProteinSequence, Set<String>>()

        val bamKmers = kmerCount.filter { it.value > 0 }.keys
        val uniqueBamKmers = bamKmers intersect hlaKmers.uniqueKmers()


        logger.info("Candidate\tMinCoverage\tMaxCoverage\tTotalCoverage")
        for ((proteinSequence, kmers) in hlaKmers.sequences()) {
            val uniqueKmers = kmers intersect uniqueBamKmers

            if (uniqueKmers.isNotEmpty() && bamKmers.containsAll(uniqueKmers)) {
                result[proteinSequence] = kmers

                val totalMatches = kmerCount.filter { it.key in kmers }.map { it.value }.sum()
                val minMatch = kmerCount.filter { it.key in kmers }.map { it.value }.min()
                val maxMatch = kmerCount.filter { it.key in kmers }.map { it.value }.max()
                logger.info("${proteinSequence.allele.toString().padEnd(14, ' ')}\t$minMatch\t$maxMatch\t$totalMatches")
                for (kmer in uniqueKmers) {
                    logger.info("       Kmer $kmer count -> ${kmerCount[kmer]}")
                }
//
//                for (kmer in uniqueKmers) {
//                    logger.info("   includes unique kmer $kmer -> ${kmerCount[kmer]}")
//                }
            } else {
//                logger.info("${proteinSequence.allele.toString().padEnd(14, ' ')} has no unique kmers")
            }
        }

        return result
    }

    private fun calculateCoverage(minMatch: Int, bamFile: String, coverage: Collection<ProteinCoverage>) {
        val bamCoverage = BamCoverage(minMatch, coverage)
        val futures = mutableListOf<Future<*>>()
        SamReaderFactory.makeDefault().open(File(bamFile)).use {

            val sequenceIndex = it.fileHeader.sequenceDictionary.getSequence("6").sequenceIndex
            val hlaa = QueryInterval(sequenceIndex, 29910331, 29913232)
            val hlac = QueryInterval(sequenceIndex, 31235946, 31240848)
            val hlab = QueryInterval(sequenceIndex, 31321260, 31325935)

            logger.info("Reading BAM file $bamFile")
            for (samRecord in it.queryOverlapping(arrayOf(hlaa, hlac, hlab))) {
                futures.add(executorService.submit { bamCoverage.accept(samRecord) })
            }

        }

        logger.info("Calculating coverage")
        futures.forEach { it.get() }

    }

    private fun calculateCoverage3(minMatch: Int, bamFile: String, coverage: Collection<ProteinCoverage>) {
        logger.info("Stuff")

        val bamCoverage = BamCoverage3(minMatch, coverage.flatMap { it.coverages })
        val futures = mutableListOf<Future<*>>()
        SamReaderFactory.makeDefault().open(File(bamFile)).use {

            val sequenceIndex = it.fileHeader.sequenceDictionary.getSequence("6").sequenceIndex
            val hlaa = QueryInterval(sequenceIndex, 29910331, 29913232)
            val hlac = QueryInterval(sequenceIndex, 31235946, 31240848)
            val hlab = QueryInterval(sequenceIndex, 31321260, 31325935)

            logger.info("Reading BAM file $bamFile")
            for (samRecord in it.queryOverlapping(arrayOf(hlaa, hlac, hlab))) {
                futures.add(executorService.submit { bamCoverage.accept(samRecord) })
            }

        }

        logger.info("Calculating coverage")
        futures.forEach { it.get() }

    }

    private fun kmerCounts(bamFile: String, kmers: Set<String>): Map<String, Int> {
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
        return kmerCount.kmerCount()
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