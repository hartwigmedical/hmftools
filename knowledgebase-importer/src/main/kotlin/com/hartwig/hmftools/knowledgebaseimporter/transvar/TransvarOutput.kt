package com.hartwig.hmftools.knowledgebaseimporter.transvar

data class TransvarOutput(val gene: String, val coordinates: String, val region: String, val info: String) {
    companion object Factory {
        operator fun invoke(transvarOutput: String): TransvarOutput {
            val sections = transvarOutput.split("\t")
            return TransvarOutput(sections[2], sections[4], sections[5], sections[6])
        }
    }
}
