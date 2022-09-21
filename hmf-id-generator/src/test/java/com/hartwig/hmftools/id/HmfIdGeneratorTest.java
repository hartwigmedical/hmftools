package com.hartwig.hmftools.id;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.common.amber.ImmutableAmberPatient;

import org.junit.Test;

public class HmfIdGeneratorTest
{
    private AmberPatient mSample1;
    private AmberPatient mSample2;
    private AmberPatient mSample3;
    private AmberPatient mSample4;
    private AmberPatient mSample5;
    private AmberPatient mSample6;

    private List<AmberPatient> mInitialSamples;
    private final List<SampleData> mInitialMappings;

    public HmfIdGeneratorTest()
    {
        mSample1 = createPatient(1, "sample1");
        mSample2 = createPatient(2, "sample2");
        mSample3 = createPatient(1, "sample3");
        mSample4 = createPatient(3, "sample4");
        mSample5 = createPatient(4, "sample5");
        mSample6 = createPatient(2, "sample6");
        mInitialSamples = Lists.newArrayList(mSample1, mSample2, mSample3, mSample4);

        mInitialMappings = HmfIdGenerator.anonymize(mInitialSamples, Lists.newArrayList());
    }

    @Test
    public void testInitialSetup()
    {
        assertMapping("sample1", "HMF000001A", mInitialMappings.get(0));
        assertMapping("sample3", "HMF000001B", mInitialMappings.get(1));
        assertMapping("sample2", "HMF000002A", mInitialMappings.get(2));
        assertMapping("sample4", "HMF000003A", mInitialMappings.get(3));
    }

    private void assertMapping(final String sampleId, final String hmfSampleId, final SampleData sample)
    {
        assertEquals(sampleId, sample.SampleId);
        assertEquals(hmfSampleId, sample.hmfSampleId());
    }

    @Test
    public void testAddNewPatientToExisting()
    {
        List<AmberPatient> patients = Lists.newArrayList(mInitialSamples);
        patients.add(mSample5);

        List<SampleData> mappings = HmfIdGenerator.anonymize(patients, mInitialMappings);
        assertMapping("sample1", "HMF000001A", mappings.get(0));
        assertMapping("sample3", "HMF000001B", mappings.get(1));
        assertMapping("sample2", "HMF000002A", mappings.get(2));
        assertMapping("sample4", "HMF000003A", mappings.get(3));
        assertMapping("sample5", "HMF000004A", mappings.get(4));
    }

    @Test
    public void testAddExistingPatientToExisting()
    {
        List<AmberPatient> patients = Lists.newArrayList(mInitialSamples);
        patients.add(mSample6);

        List<SampleData> mappings = HmfIdGenerator.anonymize(patients, mInitialMappings);
        assertMapping("sample1", "HMF000001A", mappings.get(0));
        assertMapping("sample3", "HMF000001B", mappings.get(1));
        assertMapping("sample2", "HMF000002A", mappings.get(2));
        assertMapping("sample6", "HMF000002B", mappings.get(3));
        assertMapping("sample4", "HMF000003A", mappings.get(4));
    }

    @Test
    public void testKeepSamplesEvenIfNotInAmberPatients()
    {
        List<AmberPatient> patients = Lists.newArrayList(mSample1);

        List<SampleData> mappings = HmfIdGenerator.anonymize(patients, mInitialMappings);
        assertMapping("sample1", "HMF000001A", mappings.get(0));
        assertMapping("sample3", "HMF000001B", mappings.get(1));
        assertMapping("sample2", "HMF000002A", mappings.get(2));
        assertMapping("sample4", "HMF000003A", mappings.get(3));
    }

    private AmberPatient createPatient(final int patientId, final String sampleId)
    {
        return ImmutableAmberPatient.builder().patientId(patientId).sample(sampleId).build();
    }

}