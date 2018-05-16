package com.hartwig.hmftools.actionabilityAnalyzer

data class ActionabilityOutput(val sampleId: String, val actionableTreatment: ActionableTreatment) {
    companion object {
        val header = listOf("sampleId") + ActionableTreatment.header
    }

    val record: List<String> = listOf(sampleId) + actionableTreatment.record
}
