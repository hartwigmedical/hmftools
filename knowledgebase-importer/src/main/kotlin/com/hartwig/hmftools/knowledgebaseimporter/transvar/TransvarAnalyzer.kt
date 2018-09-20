package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.HgvsAnnotation
import kotlinx.coroutines.experimental.async
import kotlinx.coroutines.experimental.runBlocking
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.FileOutputStream
import java.io.InputStreamReader
import java.io.PrintStream

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
        val args = arrayOf("--noheader", "--ensembl", "--max-candidates", "1000", "-g", "1", "-m", "2", "-l")
        logger.info("Creating transvar analyzer with command: $transvarLocation $analysisMode ${args.joinToString(" ")}")
        val process = ProcessBuilder(transvarLocation, analysisMode, *args, variantFile.absolutePath).start()
        return runBlocking {
            val errorHandler = readErrorStreamAsync(process)
            val results = readOutputStreamAsync(process)
            errorHandler.await()
            results.await()
        }
    }

    private fun readErrorStreamAsync(process: Process) = async {
        InputStreamReader(process.errorStream).buffered().lineSequence().filterNot { it.isBlank() }
                .forEach { logger.warn(" Transvar error stream produced: $it") }
    }

    private fun readOutputStreamAsync(process: Process) = async {
        InputStreamReader(process.inputStream).buffered().lineSequence().map { TransvarOutput(it) }.toList()
    }
}
