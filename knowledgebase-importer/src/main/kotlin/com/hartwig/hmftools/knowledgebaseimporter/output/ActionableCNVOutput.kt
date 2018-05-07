package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableCNVOutput(private val gene: String, private val cnvType: String, private val actionability: Actionability) {
    companion object {
        val header = listOf("gene", "cnvType") + Actionability.header
    }

    val record: List<String> = listOf(gene, cnvType) + actionability.record
}
