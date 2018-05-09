package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableFusionOutput(override val event: FusionEvent, override val actionability: Actionability) :
        ActionableItem<FusionEvent> {
    companion object {
        val fusionPairHeader = listOf("fiveGene", "threeGene") + Actionability.header
        val promiscuousGeneHeader = listOf("gene") + Actionability.header
    }

    val record: List<String> = event.record + actionability.record

    override fun toString(): String {
        return event.toString()
    }
}
