package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.ImmutableAmberPatient
import org.junit.Assert.assertEquals
import org.junit.Test

class PatientAnonymizerTest {

    private val initialPassword = "0000"
    private val sample1 = createPatient(1, "sample1")
    private val sample2 = createPatient(2, "sample2")
    private val sample3 = createPatient(1, "sample3")
    private val sample4 = createPatient(3, "sample4")
    private val initialSamples = listOf(sample1, sample2, sample3, sample4)
    private val initialHashes = PatientAnonymizer(initialPassword, initialPassword).anonymize(initialSamples, listOf())

    private val sample5 = createPatient(4, "sample5")
    private val sample6 = createPatient(2, "sample6")

    @Test
    fun testInitialSetup() {
        val mapping = AnonymizedRecord(initialPassword, initialHashes, initialSamples.map { it.sample() })
        assertEquals(AnonymizedRecord("sample1", "HMF000001A"), mapping[0])
        assertEquals(AnonymizedRecord("sample2", "HMF000002A"), mapping[1])
        assertEquals(AnonymizedRecord("sample3", "HMF000001B"), mapping[2])
        assertEquals(AnonymizedRecord("sample4", "HMF000003A"), mapping[3])
    }

    @Test
    fun testAddNewPatientToExisting() {
        val samples = initialSamples + sample5
        val hashes = PatientAnonymizer(initialPassword, initialPassword).anonymize(samples, initialHashes)
        val mapping = AnonymizedRecord(initialPassword, hashes, samples.map { it.sample() })
        assertEquals(AnonymizedRecord("sample1", "HMF000001A"), mapping[0])
        assertEquals(AnonymizedRecord("sample2", "HMF000002A"), mapping[1])
        assertEquals(AnonymizedRecord("sample3", "HMF000001B"), mapping[2])
        assertEquals(AnonymizedRecord("sample4", "HMF000003A"), mapping[3])
        assertEquals(AnonymizedRecord("sample5", "HMF000004A"), mapping[4])
    }

    @Test
    fun testAddExistingPatientToExisting() {
        val samples = initialSamples + sample6
        val hashes = PatientAnonymizer(initialPassword, initialPassword).anonymize(samples, initialHashes)
        val mapping = AnonymizedRecord(initialPassword, hashes, samples.map { it.sample() })
        assertEquals(AnonymizedRecord("sample1", "HMF000001A"), mapping[0])
        assertEquals(AnonymizedRecord("sample2", "HMF000002A"), mapping[1])
        assertEquals(AnonymizedRecord("sample3", "HMF000001B"), mapping[2])
        assertEquals(AnonymizedRecord("sample4", "HMF000003A"), mapping[3])
        assertEquals(AnonymizedRecord("sample6", "HMF000002B"), mapping[4])
    }

    @Test
    fun deletedSampleIdIsStillReserved() {
        // Remove sample 4
        var samples = listOf(sample1, sample2,  sample3)
        var hashes = PatientAnonymizer(initialPassword, initialPassword).anonymize(samples, initialHashes)

        // Adds sample 5
        samples = listOf(sample1, sample2,  sample3, sample5)
        hashes = PatientAnonymizer(initialPassword, initialPassword).anonymize(samples, hashes)

        val mapping = AnonymizedRecord(initialPassword, hashes, samples.map { it.sample() })
        assertEquals(AnonymizedRecord("sample1", "HMF000001A"), mapping[0])
        assertEquals(AnonymizedRecord("sample2", "HMF000002A"), mapping[1])
        assertEquals(AnonymizedRecord("sample3", "HMF000001B"), mapping[2])
        assertEquals(AnonymizedRecord("sample5", "HMF000004A"), mapping[3])
    }

    @Test
    fun testPasswordChange() {
        val newPassword = "1234"
        val hashes = PatientAnonymizer(initialPassword, newPassword).anonymize(initialSamples, initialHashes)

        val mapping = AnonymizedRecord(newPassword, hashes, initialSamples.map { it.sample() })
        assertEquals(AnonymizedRecord("sample1", "HMF000001A"), mapping[0])
        assertEquals(AnonymizedRecord("sample2", "HMF000002A"), mapping[1])
        assertEquals(AnonymizedRecord("sample3", "HMF000001B"), mapping[2])
    }

    @Test
    fun testAmberPatientIdChanges() {
        sample1.patientId()
        val sample1 = createPatient(110, "sample1")
        val sample2 = createPatient(210, "sample2")
        val sample3 = createPatient(110, "sample3")
        val samples = listOf(sample1, sample2, sample3)
        val hashes = PatientAnonymizer(initialPassword, initialPassword).anonymize(samples, initialHashes)
        val mapping = AnonymizedRecord(initialPassword, hashes, samples.map { it.sample() })
        assertEquals(AnonymizedRecord("sample1", "HMF000001A"), mapping[0])
        assertEquals(AnonymizedRecord("sample2", "HMF000002A"), mapping[1])
        assertEquals(AnonymizedRecord("sample3", "HMF000001B"), mapping[2])
    }

    private fun createPatient(id: Int, sample: String) = ImmutableAmberPatient.builder().patientId(id).sample(sample).build()

}