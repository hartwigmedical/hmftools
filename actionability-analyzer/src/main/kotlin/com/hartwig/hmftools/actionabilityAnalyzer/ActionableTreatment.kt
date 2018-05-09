package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableItem

data class ActionableTreatment(val treatment: String, val significance: String, val level: String,
                               val evidence: List<ActionabilityEvidence>) {
    companion object {
        val header = listOf("treatment", "significance", "level", "detailed evidence")

        operator fun invoke(items: List<ActionableItem<*>>): List<ActionableTreatment> {
            return items.groupBy { it.actionability.drug }.map { (drug, variants) ->
                ActionableTreatment(drug, significance(variants), level(variants), variants.map { ActionabilityEvidence(it) })
            }
        }

        private fun significance(items: List<ActionableItem<*>>): String {
            return "todo"
        }

        private fun level(items: List<ActionableItem<*>>): String {
            return "todo"
        }
    }

    val record: List<String> = listOf(treatment, significance, level, evidence.toString())
}
