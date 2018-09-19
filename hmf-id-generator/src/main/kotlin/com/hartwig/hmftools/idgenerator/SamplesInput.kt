package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.ids.PatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId
import org.apache.logging.log4j.LogManager

private val logger = LogManager.getLogger("SamplesInput")

data class SamplesInput(val samples: List<SampleId>, val renamedPatients: Map<PatientId, PatientId>) {
    init {
        validateRenames()
    }

    /**
     * returns all samples for this patient, accounting for potential renames
     */
    fun sampleIds(patient: PatientId): Set<SampleId> {
        return patientIds(patient).flatMap { patientId -> samples.filter { it.patientId == patientId } }.toSet()
    }

    /**
     * returns all ids for this patient, accounting for potential renames
     */
    fun patientIds(patient: PatientId): Set<PatientId> {
        val alternateIds = renamedPatients.filterValues { it == canonicalId(patient) }.flatMap { it.toPair().toList() }
        return (alternateIds + patient).toSet()
    }

    fun canonicalId(patient: PatientId) = renamedPatients[patient] ?: patient

    private fun validateRenames() {
        val chainedRenames = renamedPatients.count { (patientId, canonicalId) ->
            val canonicalIsRenamed = renamedPatients.containsKey(canonicalId)
            if (canonicalIsRenamed) {
                logger.error("Canonical id ${canonicalId.id} for ${patientId.id} is also renamed to ${renamedPatients[canonicalId]!!.id}")
            }
            canonicalIsRenamed
        }
        if (chainedRenames > 0) {
            error("Renames of canonical patientIds are not allowed")
        }
    }
}
