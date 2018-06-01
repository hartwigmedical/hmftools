package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import com.hartwig.hmftools.knowledgebaseimporter.output.*

data class ActionableKbItem(val event: ActionableEvent, val actionability: Actionability) {
    fun toActionableOutput(): ActionableItem<*>? {
        return when (event) {
            is CnvEvent            -> ActionableCNVOutput(event, actionability)
            is SomaticVariantEvent -> ActionableVariantOutput(event, actionability)
            is FusionPair          -> ActionableFusionPairOutput(event, actionability)
            is PromiscuousGene     -> ActionablePromiscuousGeneOutput(event, actionability)
        }
    }
}
