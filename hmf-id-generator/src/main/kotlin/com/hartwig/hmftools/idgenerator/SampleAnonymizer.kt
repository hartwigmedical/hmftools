package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HashId
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId
import com.hartwig.hmftools.idgenerator.ids.PatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId

class SampleAnonymizer(val password: String) {
    private val generator = IdGenerator(password)
    private val patientAnonymizer = PatientAnonymizer(password)

    fun updateSampleIds(newPassword: String, samplesInput: SamplesInput,
                        previouslyAnonymizedSamples: Collection<HmfSampleId>): AnonymizedSamples {
        val anonymizedSamples = AnonymizedSamples(password, previouslyAnonymizedSamples, samplesInput)
        val updatedPatients = patientAnonymizer.anonymize(newPassword, samplesInput, anonymizedPatients(anonymizedSamples))
        val updatedCanonicalSamples = updatedCanonicalSamples(newPassword, samplesInput, anonymizedSamples, updatedPatients)
        val samplesWithUpdatedHashes = updateHashes(newPassword, samplesInput, anonymizedSamples)
        val updatedSamples = (updatedCanonicalSamples + samplesWithUpdatedHashes).toSet()
        return AnonymizedSamples(newPassword, updatedSamples, samplesInput)
    }

    private fun updatedCanonicalSamples(newPassword: String, samplesInput: SamplesInput, anonymizedSamples: AnonymizedSamples,
                                        updatedPatients: AnonymizedPatients): List<HmfSampleId> {
        return samplesInput.canonicalPatients
                .associateBy({ it }, { updateCanonicalSamples(newPassword, it, samplesInput, anonymizedSamples) })
                .mapKeys { updatedPatients[it.key]!! }
                .mapValues { (hmfPatientId, hashedSamples) -> hashedSamples.map { HmfSampleId(it, hmfPatientId) } }
                .values.flatten()
    }

    private fun updateHashes(newPassword: String, samplesInput: SamplesInput, anonymizedSamples: AnonymizedSamples): List<HmfSampleId> {
        val hashMapping = samplesInput.hashMapping(generator, IdGenerator(newPassword))
        return anonymizedSamples.map { updateHmfSampleIdHashes(it, hashMapping) }
    }

    private fun updateHmfSampleIdHashes(hmfSampleId: HmfSampleId, hashMapping: Map<Hash, Hash>): HmfSampleId {
        val newPatientHash = hashMapping[hmfSampleId.hmfPatientId.hash] ?: hmfSampleId.hmfPatientId.hash
        val newSampleHash = hashMapping[hmfSampleId.hash] ?: hmfSampleId.hash
        return hmfSampleId.updateSampleHash(newSampleHash).updatePatientHash(newPatientHash)
    }

    private fun updateCanonicalSamples(newPassword: String, patient: PatientId, samplesInput: SamplesInput,
                                       anonymizedSamples: AnonymizedSamples): Collection<HashId> {
        val previouslyAnonymizedSamples = anonymizedSamples[samplesInput.canonicalId(patient)]
        val patientSamples = relevantPatientSamples(patient, samplesInput, previouslyAnonymizedSamples)
        return generator.update(newPassword, patientSamples.map { it.id },
                                previouslyAnonymizedSamples.map { it.hashId })
    }

    //MIVO: find all sampleIds from input that are relevant for this patient
    private fun relevantPatientSamples(patient: PatientId, samplesInput: SamplesInput,
                                       anonymizedSamples: List<HmfSampleId>): Collection<SampleId> {
        val patientSamples = samplesInput.sampleIds(patient)
        val hashes = anonymizedSamples.map { it.hash }.toSet()
        return samplesInput.samples.filter { hashes.contains(generator.hash(it.id)) } + patientSamples
    }

    private fun anonymizedPatients(anonymizedSamples: Collection<HmfSampleId>) = anonymizedSamples.map { it.hmfPatientId }.toSet()
}
