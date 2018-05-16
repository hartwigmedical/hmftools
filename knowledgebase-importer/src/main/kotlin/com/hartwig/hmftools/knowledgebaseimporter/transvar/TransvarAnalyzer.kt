package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.Annotation
import java.io.*

interface TransvarAnalyzer<in T : Annotation> {
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
        val process = ProcessBuilder(transvarLocation, analysisMode, *args, variantFile.absolutePath).start()
        val inputStream = BufferedReader(InputStreamReader(process.inputStream))
        return inputStream.lineSequence().map { TransvarOutput(it) }.toList()
    }
}
