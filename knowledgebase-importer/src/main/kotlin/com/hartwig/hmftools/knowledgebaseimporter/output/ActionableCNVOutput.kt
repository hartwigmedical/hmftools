package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableCNVOutput(override val event: CnvEvent, override val actionability: Actionability) : ActionableItem<CnvEvent> {
    companion object {
        val header = listOf("gene", "cnvType") + Actionability.header
    }

    val record: List<String> = listOf(event.gene, event.type) + actionability.record

    override fun toString(): String {
        return event.toString()
    }
}
