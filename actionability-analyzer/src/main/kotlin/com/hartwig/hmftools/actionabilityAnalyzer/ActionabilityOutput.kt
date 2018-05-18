package com.hartwig.hmftools.actionabilityAnalyzer

data class ActionabilityOutput(val sampleId: String, val patientCancerType: String?, val treatmentType: String,
                               val actionableTreatment: ActionableTreatment) {
    companion object {
        val header = listOf("sampleId", "patientCancerType", "treatmentType") + ActionableTreatment.header
    }

    val record: List<String> = listOf(sampleId, patientCancerType ?: "NULL", treatmentType) + actionableTreatment.record
}
