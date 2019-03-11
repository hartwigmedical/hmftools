package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.LocalDate;
import java.util.Map;
import java.util.Set;

import com.google.common.io.Resources;

import org.junit.Test;

public class LimsFactoryTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims").getPath();

    @Test
    public void canBuildLimsFromTestData() throws IOException {
        assertNotNull(LimsFactory.fromLimsDirectory(LIMS_DIRECTORY));
    }

    @Test(expected = IOException.class)
    public void exceptionWhenJsonFileDoesNotExist() throws IOException {
        LimsFactory.fromLimsDirectory("Does not exist");
    }

    @Test
    public void readSamplesCorrectlyFromJsonFile() throws FileNotFoundException {
        final Map<String, LimsJsonSampleData> dataPerSample =
                LimsFactory.readLimsJsonSamples(LIMS_DIRECTORY + File.separator + "lims.json");

        assertEquals(2, dataPerSample.size());

        final String refSampleId = "SAMP01010003R";
        final LimsJsonSampleData refData = dataPerSample.get(refSampleId);
        assertEquals(refSampleId, refData.sampleId());
        assertNull(refData.hospitalPatientId());
        assertEquals("2016-01-03", refData.arrivalDateString());
        assertEquals("143", refData.dnaConcentration());
        assertEquals("2016-01-02", refData.samplingDateString());
        assertEquals("NA", refData.tumorPercentageString());
        assertEquals("NA", refData.primaryTumor());
        assertEquals("PREP013V23-QC037V20-SEQ008V25", refData.labProcedures());
        assertNull(refData.labRemarks());
        assertEquals("HMFregCPCT", refData.submissionSamples());

        final String tumorSampleId = "SAMP01010003T";
        final LimsJsonSampleData tumorData = dataPerSample.get(tumorSampleId);
        assertEquals(tumorSampleId, tumorData.sampleId());
        assertEquals("something", tumorData.hospitalPatientId());
        assertEquals("2016-02-05", tumorData.arrivalDateString());
        assertEquals("143", tumorData.dnaConcentration());
        assertEquals("2016-01-04", tumorData.samplingDateString());
        assertEquals("30", tumorData.tumorPercentageString());
        assertEquals("NA", refData.primaryTumor());
        assertEquals("N/A", tumorData.labProcedures());
        assertEquals("this is a test", tumorData.labRemarks());
        assertEquals("HMFregCPCT", tumorData.submissionSamples());
    }

    @Test
    public void readSubmissionsCorrectlyFromJsonFile() throws FileNotFoundException {
        final Map<String, LimsJsonSubmissionData> dataPerSubmission =
                LimsFactory.readLimsJsonSubmissions(LIMS_DIRECTORY + File.separator + "lims.json");

        assertEquals(1, dataPerSubmission.size());

        final String submission = "ABCDEF123";
        final LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        assertEquals(submission, submissionData.submission());
    }

    @Test
    public void readCorrectlyFromPreLIMSArrivalDateFile() throws IOException {
        final Map<String, LocalDate> preLIMSArrivalDates =
                LimsFactory.readPreLIMSArrivalDateCsv(LIMS_DIRECTORY + File.separator + "pre_lims_arrival_dates.csv");

        assertNull(preLIMSArrivalDates.get("SAMP01010001T"));
        assertEquals(LimsTestUtil.toDate("2017-01-01"), preLIMSArrivalDates.get("SAMP01010003R"));
        assertEquals(LimsTestUtil.toDate("2017-01-05"), preLIMSArrivalDates.get("SAMP01010004T"));
        assertNull(preLIMSArrivalDates.get("SAMP01010005T"));
        assertNull(preLIMSArrivalDates.get("SAMP01010006R"));
        assertNull(preLIMSArrivalDates.get("DoesNotExist"));
    }

    @Test
    public void readCorrectlyFromSamplesWithoutSamplingDateFile() throws IOException {
        Set<String> samplesWithoutSamplingDate =
                LimsFactory.readSamplesWithoutSamplingDateCsv(LIMS_DIRECTORY + File.separator + "samples_without_sampling_date.csv");

        assertEquals(2, samplesWithoutSamplingDate.size());
        assertTrue(samplesWithoutSamplingDate.contains("CPCT02990001T"));
        assertFalse(samplesWithoutSamplingDate.contains("Does not exist"));
    }

    @Test
    public void readCorrectlyShallowSeqPurity() throws IOException {
        Map<String, LimsShallowSeqData> shallowSeqPuritySample =
                LimsFactory.readLimsShallowSeq(LIMS_DIRECTORY + File.separator + "shallow_seq_purity.csv");
        assertEquals(3, shallowSeqPuritySample.size());

        final String sample = "CPCT02990001T";
        final LimsShallowSeqData limsShallowSeqData = shallowSeqPuritySample.get(sample);

        assertEquals(sample, limsShallowSeqData.sampleId());
        assertEquals("0.19", limsShallowSeqData.purityShallowSeq());
    }
}