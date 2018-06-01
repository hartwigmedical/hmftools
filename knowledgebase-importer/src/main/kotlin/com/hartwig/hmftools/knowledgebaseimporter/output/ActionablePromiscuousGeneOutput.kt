package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionablePromiscuousGeneOutput(override val event: PromiscuousGene, override val actionability: Actionability) :
        ActionableItem<PromiscuousGene> {
    companion object {
        val header = PromiscuousGene.header + Actionability.header
    }

    val record: List<String> = event.record + actionability.record

    override fun toString(): String {
        return event.toString()
    }
}
