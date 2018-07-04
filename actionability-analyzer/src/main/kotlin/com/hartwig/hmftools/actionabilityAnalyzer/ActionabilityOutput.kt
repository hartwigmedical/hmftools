package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.csv.CsvData

data class ActionabilityOutput(val sampleId: String, val sampleEvent: String, val gene: String, val partnerGene: String,
                               val eventType: String, val patientCancerType: String?, val treatmentType: String, val event: String,
                               val source: String, val references: String, val drug: String, val drugType: String, val cancerType: String,
                               val level: String, val hmfLevel: String, val evidenceType: String, val significance: String,
                               val hmfResponse: String, val matchRule: String) : CsvData {
    companion object {
        operator fun invoke(sampleId: String, sampleEvent: String, gene: String, partnerGene: String,
                            eventType: String, patientCancerType: String?, treatmentType: String,
                            actionableTreatment: ActionableTreatment, matchRule: String): ActionabilityOutput {
            return ActionabilityOutput(sampleId, sampleEvent, gene, partnerGene, eventType, patientCancerType, treatmentType,
                                       actionableTreatment.event, actionableTreatment.actionability.source,
                                       actionableTreatment.actionability.reference, actionableTreatment.actionability.drug.name,
                                       actionableTreatment.actionability.drug.type, actionableTreatment.actionability.cancerType,
                                       actionableTreatment.actionability.level, actionableTreatment.actionability.hmfLevel,
                                       actionableTreatment.actionability.evidenceType, actionableTreatment.actionability.significance,
                                       actionableTreatment.actionability.hmfResponse, matchRule)

        }
    }
}
