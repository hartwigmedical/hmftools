package com.hartwig.hmftools.lilac

import com.hartwig.hmftools.lilac.hla.HlaReferenceSequence
import com.hartwig.hmftools.lilac.prot.ProteinSequenceBoundaries
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

        val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
        val sequenceAFile = "${resourcesDir}/A_prot.txt"
        val sequenceBFile = "${resourcesDir}/B_prot.txt"
        val sequenceCFile = "${resourcesDir}/C_prot.txt"

//        writeUnwrappedFiles(sequenceAFile, sequenceBFile, sequenceCFile)


        val aSequences = ProteinSequenceFile.readWrappedFile("${resourcesDir}/A_prot.txt");
        val bSequences = ProteinSequenceFile.readWrappedFile("${resourcesDir}/B_prot.txt");
        val cSequences = ProteinSequenceFile.readWrappedFile("${resourcesDir}/C_prot.txt");



    }

    private fun extractKmers() {

    }

    private fun writeUnwrappedFiles(aFile: String,  bFile: String, cFile: String) {
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.aBoundaries(), aFile, aFile.replace("txt", "unwrapped.txt"))
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.bBoundaries(), bFile, bFile.replace("txt", "unwrapped.txt"))
        ProteinSequenceFile.writeUnwrappedFile(ProteinSequenceBoundaries.cBoundaries(), cFile, cFile.replace("txt", "unwrapped.txt"))
    }


    override fun close() {
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

}