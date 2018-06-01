package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableFusionPairOutput(override val event: FusionPair, override val actionability: Actionability) :
        ActionableItem<FusionPair> {
    companion object {
        val header = FusionPair.header + Actionability.header
    }

    val record: List<String> = event.record + actionability.record

    override fun toString(): String {
        return event.toString()
    }
}
