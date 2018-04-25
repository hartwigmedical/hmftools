package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableCNVOutput(private val gene: String, private val cnvType: String, private val actionability: Actionability) {
    companion object {
        val header = listOf("gene", "cnvType", "source", "drug", "cancerType", "level", "type", "significance")
    }

    val record: List<String> = listOf(gene,
                                      cnvType,
                                      actionability.source,
                                      actionability.drug,
                                      actionability.cancerType,
                                      actionability.level,
                                      actionability.type,
                                      actionability.significance)
}
