package com.hartwig.hmftools.idgenerator

import com.google.common.annotations.VisibleForTesting
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfPatientId

class PatientAnonymizer(val password: String) {
    private val generator = IdGenerator(password)

    @VisibleForTesting
    fun anonymize(newPassword: String, samplesInput: SamplesInput, anonymizedPatients: Collection<HmfPatientId>): AnonymizedPatients {
        val canonicalPatientStrings = samplesInput.canonicalPatients.map { it.id }
        val anonymizedHashes = anonymizedPatients.map { it.hash }.toSet()
        val anonymizedNonCanonicalPatients = samplesInput.nonCanonicalPatients.filter { anonymizedHashes.contains(generator.hash(it.id)) }
        val patientIdStrings = (anonymizedNonCanonicalPatients.map { it.id } + canonicalPatientStrings).toSet()
        val updatedPatients = generator.update(newPassword, patientIdStrings, anonymizedPatients.map { it.hashId }).map { HmfPatientId(it) }
        return AnonymizedPatients(newPassword, updatedPatients.toSet(), samplesInput)
    }
}
