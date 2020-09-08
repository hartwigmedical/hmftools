package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.AmberPatient
import com.hartwig.hmftools.idgenerator.anonymizedIds.AnonymizedMap
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleIdSimple

class AmberPatientAnonymizer(private val oldPassword: String, private val newPassword: String) {

    fun anonymize(amberPatients: List<AmberPatient>, existingSamples: List<HmfSampleIdSimple>): List<HmfSampleIdSimple> {
        val result = mutableListOf<HmfSampleIdSimple>()
        var maxPatientId = existingSamples.map { x -> x.patientId }.max() ?: -1

        val amberPatientIds = amberPatients.map { x -> x.patientId() }.toSet()
        for (amberPatientId in amberPatientIds) {
            val amberSamples = amberPatients.filter { it.patientId() == amberPatientId }.map { it.sample() }
            val sampleHashMap = AnonymizedMap.create(oldPassword, newPassword, amberSamples)

            val patientExistingSamples = existingSamples
                    .filter { it.hash in sampleHashMap.oldHashes() }
                    .map { it.updateHash(sampleHashMap.fromOldHash(it.hash)) }
            result.addAll(patientExistingSamples)

            val patientId = (patientExistingSamples.map { x -> x.patientId }.max() ?: ++maxPatientId)
            var maxSampleIdForPatient = patientExistingSamples.map { x -> x.sampleId }.max() ?: 0

            for (sample in amberSamples) {
                val hash = sampleHashMap.fromSample(sample)
                if (patientExistingSamples.none { x -> x.hash == hash }) {
                    result.add(HmfSampleIdSimple(patientId, ++maxSampleIdForPatient, hash))
                }
            }
        }

        return result.sorted()
    }
}