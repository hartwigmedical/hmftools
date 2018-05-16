package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableItem

data class ActionabilityEvidence(val event: String, val source: String, val level: String, val significance: String) {
    companion object {
        operator fun invoke(item: ActionableItem<*>): ActionabilityEvidence {
            val actionability = item.actionability
            return ActionabilityEvidence(item.toString(), actionability.source, actionability.level, actionability.significance)
        }
    }

    override fun toString(): String {
        return "<$event, $source, $level, $significance>"
    }
}
