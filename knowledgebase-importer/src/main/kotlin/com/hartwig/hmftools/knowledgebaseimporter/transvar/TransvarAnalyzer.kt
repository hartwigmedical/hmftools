package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.HgvsAnnotation
import org.apache.logging.log4j.LogManager
import java.io.*

private val logger = LogManager.getLogger("TransvarAnalyzer")

interface TransvarAnalyzer<in T : HgvsAnnotation> {
    val analysisMode: String
    val transvarLocation: String

    private fun createVariantFile(variants: List<T>): File {
        val tempFile = File.createTempFile("transvar", "temp")
        tempFile.deleteOnExit()
        val writer = PrintStream(FileOutputStream(tempFile), true)
        variants.forEach { writer.println("${it.transcript}\t${it.alteration}") }
        writer.close()
        return tempFile
    }

    fun analyze(variants: List<T>): List<TransvarOutput> {
        val variantFile = createVariantFile(variants)
        val args = arrayOf("--noheader", "--ensembl", "--refversion", "hg19", "--max-candidates", "1000", "-g", "1", "-m", "2", "-l")
        logger.info("Creating transvar analyzer with command: $transvarLocation $analysisMode ${args.joinToString(" ")}")
        val process = ProcessBuilder(transvarLocation, analysisMode, *args, variantFile.absolutePath).start()
        InputStreamReader(process.errorStream).readLines().forEach { throw Exception("Transvar error: $it") }
        val inputStream = BufferedReader(InputStreamReader(process.inputStream))
        return inputStream.lineSequence().map { TransvarOutput(it) }.toList()
    }
}
