package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.AmberPatient

class PatientAnonymizer {

    fun anonymize(amberPatients: List<AmberPatient>, existingSamples: List<HmfSample>): List<HmfSample> {
        val newSamples = mutableListOf<HmfSample>()
        var maxPatientId = existingSamples.map { x -> x.patientId }.max() ?: 0

        val amberPatientIds = amberPatients.map { x -> x.patientId() }.toSet()
        for (amberPatientId in amberPatientIds) {
            val allPatientSamples = amberPatients.filter { it.patientId() == amberPatientId }.map { it.sample() }
            val existingPatientSamples = existingSamples.filter { x -> allPatientSamples.contains(x.sample) }

            val hmfPatientId = (existingPatientSamples.map { x -> x.patientId }.max() ?: ++maxPatientId)
            var maxSampleIdForPatient = existingPatientSamples.map { x -> x.sampleId }.max() ?: 0

            for (sample in allPatientSamples) {
                if (existingPatientSamples.none { x -> x.sample == sample }) {
                    newSamples.add(HmfSample(hmfPatientId, ++maxSampleIdForPatient, false, sample))
                }
            }
        }
        return (existingSamples + newSamples).sorted()
    }
}