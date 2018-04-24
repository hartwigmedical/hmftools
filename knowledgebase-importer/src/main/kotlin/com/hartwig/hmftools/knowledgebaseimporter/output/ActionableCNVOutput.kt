package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableCNVOutput(private val gene: String, private val type: String, private val actionability: Actionability) {
    companion object {
        val header = listOf("gene", "type", "source", "drug", "cancerType", "level", "type", "significance")
    }

    val record: List<String> = listOf(gene,
                                      type,
                                      actionability.source,
                                      actionability.drug,
                                      actionability.cancerType,
                                      actionability.level,
                                      actionability.type,
                                      actionability.significance)
}
