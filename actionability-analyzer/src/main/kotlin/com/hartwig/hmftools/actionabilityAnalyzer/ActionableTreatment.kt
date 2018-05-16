package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableItem

data class ActionableTreatment(val event: String, val actionability: Actionability) {
    companion object {
        val header = listOf("event") + Actionability.header

        operator fun invoke(items: List<ActionableItem<*>>): List<ActionableTreatment> {
            return items.map { ActionableTreatment(it.event.toString(), it.actionability) }
        }
    }

    val record: List<String> = listOf(event) + actionability.record
}
