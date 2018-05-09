package com.hartwig.hmftools.knowledgebaseimporter.output

data class CnvEvent(val gene: String, val type: String) : ActionableEvent {
    companion object {
        val header = listOf("gene", "type")
    }

    val record: List<String> = listOf(gene, type)

    override fun toString(): String {
        return "$gene $type"
    }
}
