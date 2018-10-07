package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.ids.CanonicalPatientId
import com.hartwig.hmftools.idgenerator.ids.PatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId
import org.apache.logging.log4j.LogManager

data class SamplesInput(val samples: List<SampleId>, val patientsMap: Map<PatientId, PatientId> = emptyMap()) {
    private val logger = LogManager.getLogger(this.javaClass)
    private val allPatientIdsPerPatientId: Map<PatientId, Set<PatientId>> = allPatientIdsPerPatientId()
    private val allSampleIdsPerPatientId: Map<PatientId, Set<SampleId>> = allSamplesPerPatientId()

    init {
        validatePatientsMap()
    }

    val canonicalPatients: Set<CanonicalPatientId> =
            (samples.map { canonicalId(it.patientId) } + patientsMap.values.map { CanonicalPatientId(it) }).toSet()
    val nonCanonicalPatients: Set<PatientId> = patientsMap.keys

    /**
     * returns all samples for this patient, accounting for potential renames
     */
    fun sampleIds(patient: PatientId): Set<SampleId> = allSampleIdsPerPatientId.getValue(patient)

    /**
     * returns all ids for this patient, accounting for potential renames
     */
    fun patientIds(patient: PatientId): Set<PatientId> = allPatientIdsPerPatientId.getValue(patient)

    fun canonicalId(patient: PatientId) = CanonicalPatientId(patientsMap[patient] ?: patient)

    fun hashMapping(generator: IdGenerator, newGenerator: IdGenerator): Map<Hash, Hash> {
        val samplePlaintexts = samples.map { it.id }
        val patientPlaintexts = samples.map { it.patientId.id } + patientsMap.flatMap { it.toPair().toList() }.map { it.id }
        val allPlaintexts = samplePlaintexts + patientPlaintexts
        return allPlaintexts.associateBy({ generator.hash(it) }, { newGenerator.hash(it) })
    }

    private fun allPatientIdsPerPatientId(): Map<PatientId, Set<PatientId>> {
        val allPatientIds = samples.map { it.patientId } + patientsMap.flatMap { it.toPair().toList() }
        val canonicalAlternateIds = patientsMap.map { it.toPair() }
                .groupBy({ CanonicalPatientId(it.second) }) { it.first }.withDefault { _ -> emptyList() }
        return allPatientIds.associateBy({ it }) { patientId ->
            val canonicalId = canonicalId(patientId)
            (canonicalAlternateIds.getValue(canonicalId) + canonicalId.patientId).toSet()
        }.withDefault { setOf(it) }
    }

    private fun allSamplesPerPatientId(): Map<PatientId, Set<SampleId>> {
        val samplesPerPatientId = samples.groupBy { it.patientId }.withDefault { emptyList() }
        return allPatientIdsPerPatientId.mapValues { it.value.flatMap { patientId -> samplesPerPatientId.getValue(patientId) }.toSet() }
                .withDefault { emptySet() }
    }

    private fun validatePatientsMap() {
        val chainedCanonicalMappings = patientsMap.count { (patientId, canonicalId) ->
            val canonicalIsRenamed = patientsMap.containsKey(canonicalId)
            if (canonicalIsRenamed) {
                logger.error("Canonical id ${canonicalId.id} for ${patientId.id} is also renamed to ${patientsMap[canonicalId]!!.id}")
            }
            canonicalIsRenamed
        }
        if (chainedCanonicalMappings > 0) {
            error("Renames of canonical patientIds are not allowed")
        }
    }
}
