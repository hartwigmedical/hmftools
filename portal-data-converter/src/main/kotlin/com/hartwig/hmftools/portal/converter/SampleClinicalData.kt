package com.hartwig.hmftools.portal.converter

import org.apache.commons.csv.CSVRecord
import java.time.LocalDate
import java.time.format.DateTimeFormatter

data class SampleClinicalData(val sampleId: String, val gender: String, val ageAtEnrollment: String, val cancerType: String,
                              val specimenType: String, val specimenTypeOther: String) {
    companion object Factory {
        private val DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd")

        operator fun invoke(record: CSVRecord): SampleClinicalData {
            val gender = determineGender(record.get(GENDER_FIELD))
            val ageAtEnrollment = determineAgeAtEnrollment(record.get(BIRTH_YEAR_FIELD), record.get(REGISTRATION_DATE_FIELD))
            val cancerType = determineCancerType(record.get(CANCER_TYPE_FIELD))
            val (specimenType, specimenTypeOther) = determineSpecimenType(record.get(BIOPSY_SITE_FIELD))
            return SampleClinicalData(record.get(SAMPLE_ID_FIELD), gender, ageAtEnrollment, cancerType, specimenType, specimenTypeOther)
        }

        private fun determineGender(ecrfGender: String?): String {
            return when (ecrfGender?.trim()) {
                "male" -> "1"
                "female" -> "2"
                else -> DEFAULT_VALUE
            }
        }

        private fun determineCancerType(ecrfCancerType: String?): String {
            return if (ecrfCancerType == null || ecrfCancerType.isBlank()) "Missing" else ecrfCancerType
        }

        private fun determineSpecimenType(ecrfBiopsySite: String?): Pair<String, String> {
            return when (ecrfBiopsySite?.trim()?.toLowerCase()) {
                null -> Pair(DEFAULT_VALUE, DEFAULT_VALUE)
                "primary" -> Pair("Primary tumour - solid tissue", DEFAULT_VALUE)
                "lymph node" -> Pair("Metastatic tumour - lymph node", DEFAULT_VALUE)
                else -> Pair("Metastatic tumour - other", ecrfBiopsySite)
            }
        }

        private fun determineAgeAtEnrollment(birthYearString: String?, registrationDateString: String?): String {
            return if (birthYearString == null || registrationDateString == null) {
                DEFAULT_VALUE
            } else {
                try {
                    val birthYear = birthYearString.toInt()
                    val registrationDate = LocalDate.parse(registrationDateString, DATE_FORMATTER)
                    (registrationDate.year - birthYear).toString()
                } catch (e: Exception) {
                    DEFAULT_VALUE
                }
            }
        }
    }
}

