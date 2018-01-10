package com.hartwig.hmftools.portal.converter.records

import com.hartwig.hmftools.portal.converter.SampleClinicalData
import com.hartwig.hmftools.portal.converter.records.donor.Donor
import com.hartwig.hmftools.portal.converter.records.sample.Sample
import com.hartwig.hmftools.portal.converter.records.specimen.Specimen

data class SampleRecords(private val patientId: String, private val sampleId: String, private val clinicalData: SampleClinicalData) {
    val donor = Donor(patientId, clinicalData.gender, clinicalData.ageAtEnrollment)
    val sample = Sample(sampleId)
    val specimen = Specimen(sampleId, clinicalData.specimenType, clinicalData.specimenTypeOther, patientId)
}
