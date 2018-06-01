package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.csv.CsvData

data class ActionabilityOutput(val sampleId: String, val patientCancerType: String?, val treatmentType: String,
                               val actionableTreatment: ActionableTreatment) : CsvData
