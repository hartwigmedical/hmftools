package com.hartwig.hmftools.idgenerator.ids

data class CanonicalPatientId(val patientId: PatientId) {
    val id = patientId.id
}
