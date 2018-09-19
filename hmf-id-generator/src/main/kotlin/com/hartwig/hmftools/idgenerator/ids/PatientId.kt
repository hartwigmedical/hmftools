package com.hartwig.hmftools.idgenerator.ids

import org.apache.logging.log4j.LogManager

data class PatientId(val type: IdType, val id: String) {
    companion object {
        private val logger = LogManager.getLogger(this::class.java)

        operator fun invoke(id: String): PatientId? {
            return if (id.length == 12 && IdType.values().any { id.startsWith(it.name) }) {
                PatientId(IdType.valueOf(id.substring(0, 4)), id)
            } else {
                logger.warn("Could not create patient id from: $id")
                return null
            }
        }
    }
}
