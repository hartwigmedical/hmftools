package com.hartwig.hmftools.lilac

import com.hartwig.hmftools.lilac.hla.HlaReferenceSequence
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile
import org.apache.logging.log4j.LogManager

fun main(args: Array<String>) {
    LilacApplication().use { x -> x.run() }
}


class LilacApplication : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis()

    override fun run() {

//        logger.info("Reading HLA Sequences")
        val allSequences = mutableListOf<HlaReferenceSequence>()
//        val sequenceA = readHlaSequence("/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/A_prot.fasta")
//        val sequenceB = readHlaSequence("/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/B_prot.fasta")
//        val sequenceC = readHlaSequence("/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/C_prot.fasta")


        ProteinSequenceFile.writeUnwrappedFile("/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/A_prot.txt", "/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/A_prot.unwrapped.txt")

//        allSequences.addAll(sequenceA)
//        allSequences.addAll(sequenceB)
//        allSequences.addAll(sequenceC)

//        KmerByGeneCount("/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/KmerByGeneCount.tmp.txt").doIt(sequenceA, sequenceB, sequenceC)
//        HlaUniqueKmerCountSuperSlow("/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/HlaUniqueKmerCount.tmp.txt").doIt(sequenceA, sequenceB, sequenceC)
//        HlaUniqueKmerCountFast(10, 400, "/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/HlaUniqueKmerCount.10.400.txt").doIt(sequenceA, sequenceB, sequenceC)


//        logger.info("Extracting kmers")
//        val kmers = sequences.flatMap { x -> x.rollingKmers(10) }.distinct().sorted()

//        for (i in 0..200) {
//            val jon = UniqueKmers.uniqueKmers(sequences[i], sequences)
//            println("${sequences[i].allele.fourDigitName()} -> ${jon[0].length} -> ${jon[0]}")
//        }
//
//
//
////        val jon = UniqueKmers.uniqueKmers(sequences[3], sequences)
////        println(jon[0].length)
//
////        val file = KmerCountFile("/Users/jon/hmf/repos/hmftools/lilac/src/main/resources/A_count.txt", sequences)
////        file.addEntry(kmers)
//
//
//        println(kmers.size)

    }

    private fun extractKmers() {

    }


    override fun close() {
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

}