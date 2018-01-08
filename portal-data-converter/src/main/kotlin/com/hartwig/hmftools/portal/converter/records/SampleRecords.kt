package com.hartwig.hmftools.portal.converter.records

import com.hartwig.hmftools.portal.converter.SampleClinicalData
import com.hartwig.hmftools.portal.converter.ids.extractPatientId
import com.hartwig.hmftools.portal.converter.records.donor.Donor
import com.hartwig.hmftools.portal.converter.records.sample.Sample
import com.hartwig.hmftools.portal.converter.records.specimen.Specimen

data class SampleRecords(private val clinicalData: SampleClinicalData) {
    private val patientId = extractPatientId(clinicalData.sampleId)
    val donor = Donor(patientId, clinicalData.gender, clinicalData.ageAtEnrollment)
    val sample = Sample(clinicalData.sampleId)
    val specimen = Specimen(clinicalData.sampleId, clinicalData.specimenType, clinicalData.specimenTypeOther, patientId)
}