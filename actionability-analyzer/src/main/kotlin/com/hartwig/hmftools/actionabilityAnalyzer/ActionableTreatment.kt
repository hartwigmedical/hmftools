package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableItem
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability

data class ActionableTreatment(val event: String, val actionability: Actionability) : CsvData {
    companion object {
        operator fun invoke(items: List<ActionableItem<*>>): List<ActionableTreatment> {
            val actionabilityMap = groupedReferencesActionability(items.map { it.actionability })
            return items.map { ActionableTreatment(it.event.eventString(), actionabilityMap[actionabilityKey(it.actionability)]!!) }
        }

        //MIVO: map multiple identical actionability items that only differ on reference into single item with concatenated references
        private fun groupedReferencesActionability(actionabilityItems: List<Actionability>): Map<Actionability, Actionability> {
            return actionabilityItems.groupBy { actionabilityKey(it) }.mapValues { (key, actionabilityList) ->
                val concatenatedReferences = actionabilityList.map { it.reference }.toSet().joinToString(",")
                key.copy(reference = concatenatedReferences)
            }
        }

        private fun actionabilityKey(actionability: Actionability): Actionability = actionability.copy(reference = "")
    }
}
