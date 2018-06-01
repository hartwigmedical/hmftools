package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableItem
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability

data class ActionableTreatment(val event: String, val actionability: Actionability) : CsvData {
    companion object {

        operator fun invoke(items: List<ActionableItem<*>>): List<ActionableTreatment> {
            return items.map { ActionableTreatment(it.event.eventString(), it.actionability) }
        }
    }
}
