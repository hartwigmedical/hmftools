package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarAnalyzer.Companion.commonArgs
import java.io.BufferedReader
import java.io.InputStreamReader

data class TransvarCdnaAnalyzer(private val transvarLocation: String) : TransvarAnalyzer {

    override fun analyze(variants: List<Pair<Transcript, Impact>>): List<TransvarOutput> {
        val variantFile = createVariantFile(variants)
        val process = ProcessBuilder(transvarLocation, "canno", *commonArgs, variantFile.absolutePath).start()
        val inputStream = BufferedReader(InputStreamReader(process.inputStream))
        return inputStream.lineSequence().map { TransvarOutput(it) }.toList()
    }
}
