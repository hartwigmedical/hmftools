package com.hartwig.hmftools.portal.converter

import com.hartwig.hmftools.common.ecrf.projections.PortalClinicalData
import com.hartwig.hmftools.portal.converter.records.donor.Donor
import com.hartwig.hmftools.portal.converter.records.sample.Sample
import com.hartwig.hmftools.portal.converter.records.specimen.Specimen
import java.time.LocalDate

data class SampleClinicalData(val cpctId: String, val sampleId: String, private val gender: String, private val ageAtEnrollment: String,
                              val cancerType: String, private val specimenType: String, private val specimenTypeOther: String,
                              val hmfId: String, val sampleHmfId: String) {
    companion object Factory {
        operator fun invoke(patientData: PortalClinicalData): SampleClinicalData {
            val gender = determineGender(patientData.gender())
            val ageAtEnrollment = determineAgeAtEnrollment(patientData.birthYear(), patientData.registrationDate())
            val cancerType = determineCancerType(patientData.cancerType())
            val (specimenType, specimenTypeOther) = determineSpecimenType(patientData.biopsyType())
            return SampleClinicalData(patientData.patientIdentifier(),
                                      patientData.sampleId(),
                                      gender,
                                      ageAtEnrollment,
                                      cancerType,
                                      specimenType,
                                      specimenTypeOther,
                                      patientData.hmfId(),
                                      patientData.sampleHmfId())
        }

        private fun determineGender(ecrfGender: String): String {
            return when (ecrfGender.trim()) {
                "male"   -> "1"
                "female" -> "2"
                else     -> DEFAULT_VALUE
            }
        }

        private fun determineCancerType(ecrfCancerType: String): String {
            return if (ecrfCancerType.isBlank()) "Missing" else ecrfCancerType
        }

        private fun determineSpecimenType(ecrfBiopsySite: String): Pair<String, String> {
            return when {
                ecrfBiopsySite.isBlank() -> Pair("Biopsy site: -", DEFAULT_VALUE)
                else                     -> Pair("Biopsy site: " + ecrfBiopsySite, DEFAULT_VALUE)
            }
        }

        private fun determineAgeAtEnrollment(birthYearString: String, registrationDateString: String): String {
            return if (birthYearString.isBlank() || registrationDateString.isBlank()) {
                DEFAULT_VALUE
            } else {
                try {
                    val birthYear = birthYearString.toInt()
                    val registrationDate = LocalDate.parse(registrationDateString, PortalClinicalData.FORMATTER)
                    (registrationDate.year - birthYear).toString()
                } catch (e: Exception) {
                    DEFAULT_VALUE
                }
            }
        }
    }

    val donor = Donor(hmfId, gender, ageAtEnrollment)
    val sample = Sample(sampleHmfId)
    val specimen = Specimen(sampleHmfId, specimenType, specimenTypeOther, hmfId)
}
