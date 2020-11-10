package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.ImmutableAmberPatient
import org.junit.Assert.assertEquals
import org.junit.Test

class PatientAnonymizerTest {

    private val sample1 = createPatient(1, "sample1")
    private val sample2 = createPatient(2, "sample2")
    private val sample3 = createPatient(1, "sample3")
    private val sample4 = createPatient(3, "sample4")
    private val initialSamples = listOf(sample1, sample2, sample3, sample4)
    private val initialMapping = PatientAnonymizer().anonymize(initialSamples, listOf())

    private val sample5 = createPatient(4, "sample5")
    private val sample6 = createPatient(2, "sample6")

    @Test
    fun testInitialSetup() {
        assertMapping("sample1", "HMF000001A", initialMapping[0])
        assertMapping("sample3", "HMF000001B", initialMapping[1])
        assertMapping("sample2", "HMF000002A", initialMapping[2])
        assertMapping("sample4", "HMF000003A", initialMapping[3])
    }

    private fun assertMapping(sample: String, hmfSample: String, victim: HmfSample) {
        assertEquals(sample, victim.sample)
        assertEquals(hmfSample, victim.hmfSample)

    }

    @Test
    fun testAddNewPatientToExisting() {
        val samples = initialSamples + sample5
        val mapping = PatientAnonymizer().anonymize(samples, initialMapping)
        assertMapping("sample1", "HMF000001A", mapping[0])
        assertMapping("sample3", "HMF000001B", mapping[1])
        assertMapping("sample2", "HMF000002A", mapping[2])
        assertMapping("sample4", "HMF000003A", mapping[3])
        assertMapping("sample5", "HMF000004A", mapping[4])
    }

    @Test
    fun testAddExistingPatientToExisting() {
        val samples = initialSamples + sample6
        val mapping = PatientAnonymizer().anonymize(samples, initialMapping)
        assertMapping("sample1", "HMF000001A", mapping[0])
        assertMapping("sample3", "HMF000001B", mapping[1])
        assertMapping("sample2", "HMF000002A", mapping[2])
        assertMapping("sample6", "HMF000002B", mapping[3])
        assertMapping("sample4", "HMF000003A", mapping[4])
    }

    @Test
    fun testKeepSamplesEvenIfNotInAmberPatients() {
        val samples = listOf(sample1)
        val mapping = PatientAnonymizer().anonymize(samples, initialMapping)
        assertMapping("sample1", "HMF000001A", mapping[0])
        assertMapping("sample3", "HMF000001B", mapping[1])
        assertMapping("sample2", "HMF000002A", mapping[2])
        assertMapping("sample4", "HMF000003A", mapping[3])
    }


    private fun createPatient(id: Int, sample: String) = ImmutableAmberPatient.builder().patientId(id).sample(sample).build()

}