package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableCNVOutput(override val event: CnvEvent, override val actionability: Actionability) : ActionableItem<CnvEvent> {
    companion object {
        val header = CnvEvent.header + Actionability.header
    }

    val record: List<String> = event.record + actionability.record

    override fun toString(): String {
        return event.toString()
    }
}
