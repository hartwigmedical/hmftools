package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableFusionOutput(override val event: FusionEvent, override val actionability: Actionability) :
        ActionableItem<FusionEvent> {
    companion object {
        val fusionPairHeader = FusionPair.header + Actionability.header
        val promiscuousGeneHeader = PromiscuousGene.header + Actionability.header
    }

    val record: List<String> = event.record + actionability.record

    override fun toString(): String {
        return event.toString()
    }
}
