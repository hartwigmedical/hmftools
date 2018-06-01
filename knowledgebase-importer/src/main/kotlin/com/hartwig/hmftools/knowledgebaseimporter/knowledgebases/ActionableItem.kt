package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import com.hartwig.hmftools.knowledgebaseimporter.output.*

interface ActionableItem<T : ActionableEvent> {
    val event: T
    val actionability: Actionability

    companion object {
        operator fun invoke(event: ActionableEvent, actionability: Actionability): ActionableItem<*> {
            return when (event) {
                is CnvEvent            -> ActionableCNVOutput(event, actionability)
                is SomaticVariantEvent -> ActionableVariantOutput(event, actionability)
                is FusionPair          -> ActionableFusionPairOutput(event, actionability)
                is PromiscuousGene     -> ActionablePromiscuousGeneOutput(event, actionability)
            }
        }
    }
}
