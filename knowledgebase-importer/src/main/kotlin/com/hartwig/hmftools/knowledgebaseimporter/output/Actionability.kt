package com.hartwig.hmftools.knowledgebaseimporter.output

data class Actionability(val source: String, val cancerType: String, val drug: String, val level: String, val significance: String,
                         val evidenceType: String) {
    companion object {
        val header = listOf("source", "drug", "cancerType", "level", "evidenceType", "significance")
    }

    val record: List<String> = listOf(source, drug, cancerType, level, evidenceType, significance)
}
