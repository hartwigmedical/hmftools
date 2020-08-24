package com.hartwig.hmftools.paddle.cohort

import com.hartwig.hmftools.common.amber.ImmutableAmberPatient
import com.hartwig.hmftools.common.purple.purity.ImmutableSamplePurity
import junit.framework.Assert.assertEquals
import org.junit.Test

class HighestPuritySampleTest {

    @Test
    fun testHighestPurity() {
        val patientId = 1
        val sampleId1 = "1"
        val sampleId2 = "2"

        val amberPatient1 = ImmutableAmberPatient.builder().patientId(patientId).sample(sampleId1).build();
        val amberPatient2 = ImmutableAmberPatient.builder().patientId(patientId).sample(sampleId2).build();
        val patients = mutableListOf(amberPatient1, amberPatient2)
        patients.shuffle()

        val samplePurity1 = ImmutableSamplePurity.builder().sampleId(sampleId1).purity(0.9).build();
        val samplePurity2 = ImmutableSamplePurity.builder().sampleId(sampleId2).purity(0.8).build()
        val purities = mutableListOf(samplePurity1, samplePurity2)
        purities.shuffle()

        val highestPurityCohort = HighestPuritySample.highestPurityCohort(purities, patients)
        assertEquals(1, highestPurityCohort.size)
        assertEquals(sampleId1, highestPurityCohort[0].sampleId)
        assertEquals(0.9, highestPurityCohort[0].purity)
    }

}