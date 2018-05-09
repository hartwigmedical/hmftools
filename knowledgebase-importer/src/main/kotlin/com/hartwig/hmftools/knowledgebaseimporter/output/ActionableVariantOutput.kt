package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableVariantOutput(val gene: String, override val event: SomaticVariantEvent, override val actionability: Actionability) :
        ActionableItem<SomaticVariantEvent> {
    companion object {
        val header = listOf("gene") + SomaticVariantEvent.header + Actionability.header
    }

    val record: List<String> = listOf(gene) + event.record + actionability.record

    override fun toString(): String {
        return "$gene $event"
    }
}
