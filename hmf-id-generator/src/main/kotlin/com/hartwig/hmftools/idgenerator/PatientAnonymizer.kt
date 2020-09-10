package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.AmberPatient
import com.hartwig.hmftools.idgenerator.anonymizedIds.AnonymizedMap
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId
import org.apache.logging.log4j.LogManager

private val logger = LogManager.getLogger("PatientAnonymizer")

class PatientAnonymizer(private val oldPassword: String, private val newPassword: String) {

    fun anonymize(amberPatients: List<AmberPatient>, existingSamples: List<HmfSampleId>): List<HmfSampleId> {
        val newSamples = mutableListOf<HmfSampleId>()
        val oldSamples = mutableListOf<HmfSampleId>()
        val missingSamples = mutableListOf<HmfSampleId>()
        var maxPatientId = existingSamples.map { x -> x.patientId }.max() ?: 0

        val amberPatientIds = amberPatients.map { x -> x.patientId() }.toSet()
        for (amberPatientId in amberPatientIds) {
            val amberSamples = amberPatients.filter { it.patientId() == amberPatientId }.map { it.sample() }
            val sampleHashMap = AnonymizedMap.create(oldPassword, newPassword, amberSamples)

            val patientExistingSamples = existingSamples
                    .filter { it.hash in sampleHashMap.oldHashes() }
                    .map { it.updateHash(sampleHashMap.fromOldHash(it.hash)) }
            oldSamples.addAll(patientExistingSamples)

            val patientId = (patientExistingSamples.map { x -> x.patientId }.max() ?: ++maxPatientId)
            var maxSampleIdForPatient = patientExistingSamples.map { x -> x.sampleId }.max() ?: 0

            for (sample in amberSamples) {
                val hash = sampleHashMap.fromSample(sample)
                if (patientExistingSamples.none { x -> x.hash == hash }) {
                    newSamples.add(HmfSampleId(patientId, ++maxSampleIdForPatient, hash))
                }
            }
        }

        // Check for incorrect password
        if (existingSamples.isNotEmpty() && oldSamples.isEmpty()) {
            logger.error("None of the existing hashes matched the supplied samples.")
            throw IllegalStateException("Incorrect password")
        }

        // Detect deleted samples
        val oldUpdatedSamples = oldSamples.map { x -> x.plaintext }.toSet()
        for (existingSample in existingSamples) {
            if (!oldUpdatedSamples.contains(existingSample.plaintext)) {
                missingSamples.add(existingSample.updateHash(Hash("na")))
            }
        }

        return (newSamples + oldSamples + missingSamples).sorted()
    }
}