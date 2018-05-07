package com.hartwig.hmftools.knowledgebaseimporter.output

data class Actionability(private val source: String, private val cancerType: String, private val drug: String, private val level: String,
                         private val significance: String, private val evidenceType: String) {
    companion object {
        val header = listOf("source", "drug", "cancerType", "level", "evidenceType", "significance")
    }

    val record: List<String> = listOf(source, drug, cancerType, level, evidenceType, significance)
}