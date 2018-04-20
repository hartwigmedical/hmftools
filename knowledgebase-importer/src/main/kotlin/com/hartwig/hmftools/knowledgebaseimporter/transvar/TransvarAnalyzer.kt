package com.hartwig.hmftools.knowledgebaseimporter.transvar

import java.io.File
import java.io.FileOutputStream
import java.io.PrintStream

typealias Transcript = String
typealias Impact = String

interface TransvarAnalyzer {
    companion object {
        val commonArgs = arrayOf("--noheader", "--ensembl", "--max-candidates", "1000", "-g", "1", "-m", "2", "-l")
    }

    fun createVariantFile(variants: List<Pair<Transcript, Impact>>): File {
        val tempFile = File.createTempFile("transvar", "temp")
        tempFile.deleteOnExit()
        val writer = PrintStream(FileOutputStream(tempFile), true)
        variants.forEach { writer.println("${it.first}\t${it.second}") }
        writer.close()
        return tempFile
    }

    fun analyze(variants: List<Pair<Transcript, Impact>>): List<TransvarOutput>

}
