package com.hartwig.hmftools.idgenerator.ids

import org.apache.logging.log4j.LogManager

data class SampleId(val patientId: PatientId, val id: String) {
    companion object {
        private const val SAMPLE_ID_MIN_SIZE = 13

        private val logger = LogManager.getLogger(this::class.java)

        operator fun invoke(sampleId: String): SampleId? {
            return if (sampleId.length >= SAMPLE_ID_MIN_SIZE) {

                val patientId = PatientId(extractPatientIdString(sampleId))
                if (patientId == null) null
                else SampleId(patientId, sampleId)
            } else {
                logger.warn("Could not create sample id from: $sampleId")
                return null
            }
        }

        private fun extractPatientIdString(sampleId: String): String {
            val pattern = "([A-Z]+[0-9]+)".toRegex()
            return pattern.find(sampleId)?.groupValues?.get(1) ?: ""
        }
    }
}
