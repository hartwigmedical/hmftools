package com.hartwig.hmftools.knowledgebaseimporter.output

data class CnvEvent(val gene: String, val cnvType: String) : ActionableEvent {
    companion object {
        val header = listOf("gene", "cnvType")
    }

    val record: List<String> = listOf(gene, cnvType)

    override fun toString(): String {
        return "$gene $cnvType"
    }
}
