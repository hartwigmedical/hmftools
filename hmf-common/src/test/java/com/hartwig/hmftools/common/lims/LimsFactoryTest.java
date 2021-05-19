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
    public void canCreateEmptyLims() {
        assertNotNull(LimsFactory.empty());
    }

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
        Map<String, LimsJsonSampleData> dataPerSampleBarcode =
                LimsFactory.readLimsJsonSamples(LIMS_DIRECTORY + File.separator + "lims.json");

        assertEquals(4, dataPerSampleBarcode.size());

        String refSampleBarcode = "FR1234";
        LimsJsonSampleData refData = dataPerSampleBarcode.get(refSampleBarcode);
        assertNull(refData.hospitalPatientId());
        assertEquals("SAMP01010003R", refData.sampleId());
        assertEquals("2016-01-03", refData.arrivalDate());
        assertEquals("143", refData.dnaConcentration());
        assertEquals("2016-01-02", refData.samplingDate());
        assertEquals("NA", refData.pathologyTumorPercentage());
        assertEquals("NA", refData.primaryTumor());
        assertEquals("PREP013V23-QC037V20-SEQ008V25", refData.labProcedures());
        assertEquals("HMFregCPCT", refData.submission());

        String tumorSampleBarcode = "FR1236";
        LimsJsonSampleData tumorData = dataPerSampleBarcode.get(tumorSampleBarcode);
        assertEquals(tumorSampleBarcode, tumorData.tumorBarcode());
        assertEquals("SAMP01010003T", tumorData.sampleId());
        assertEquals("something", tumorData.hospitalPatientId());
        assertEquals("2016-02-05", tumorData.arrivalDate());
        assertEquals("143", tumorData.dnaConcentration());
        assertEquals("2016-01-04", tumorData.samplingDate());
        assertEquals("30", tumorData.pathologyTumorPercentage());
        assertEquals("NA", tumorData.primaryTumor());
        assertEquals("N/A", tumorData.labProcedures());
        assertEquals("PREPV-QC037V20-SEQ008V25", tumorData.labSopVersions());
        assertEquals("HMFregCPCT", tumorData.submission());
        assertEquals("", tumorData.hospitalPathologySampleId());
        assertEquals("", tumorData.germlineReportingLevel());
    }

    @Test
    public void readSubmissionsCorrectlyFromJsonFile() throws FileNotFoundException {
        Map<String, LimsJsonSubmissionData> dataPerSubmission =
                LimsFactory.readLimsJsonSubmissions(LIMS_DIRECTORY + File.separator + "lims.json");

        assertEquals(1, dataPerSubmission.size());

        String submission = "ABCDEF123";
        LimsJsonSubmissionData submissionData = dataPerSubmission.get(submission);
        assertEquals(submission, submissionData.submission());
    }

    @Test
    public void readCorrectlyFromPreLimsArrivalDateFile() throws IOException {
        Map<String, LocalDate> preLimsArrivalDates =
                LimsFactory.readPreLimsArrivalDateTsv(LIMS_DIRECTORY + File.separator + "pre_lims_arrival_dates.tsv");

        assertNull(preLimsArrivalDates.get("SAMP01010001T"));
        assertEquals(LimsTestUtil.toDate("2017-01-01"), preLimsArrivalDates.get("SAMP01010003R"));
        assertEquals(LimsTestUtil.toDate("2017-01-05"), preLimsArrivalDates.get("SAMP01010004T"));
        assertNull(preLimsArrivalDates.get("SAMP01010005T"));
        assertNull(preLimsArrivalDates.get("SAMP01010006R"));
        assertNull(preLimsArrivalDates.get("DoesNotExist"));
    }

    @Test
    public void readCorrectlyFromSamplesWithoutSamplingDateFile() throws IOException {
        Set<String> samplesWithoutSamplingDate =
                LimsFactory.readSingleColumnTsv(LIMS_DIRECTORY + File.separator + "samples_without_sampling_date.tsv");

        assertEquals(2, samplesWithoutSamplingDate.size());
        assertTrue(samplesWithoutSamplingDate.contains("SAMP01011234T"));
        assertFalse(samplesWithoutSamplingDate.contains("Does not exist"));
    }

    @Test
    public void readCorrectlyFromShallowSeqPurityFile() throws IOException {
        Map<String, LimsShallowSeqData> shallowSeqPuritySample =
                LimsFactory.readLimsShallowSeqTsv(LIMS_DIRECTORY + File.separator + "shallow_seq_purity.tsv");
        assertEquals(3, shallowSeqPuritySample.size());

        String barcode = "FR1";

        LimsShallowSeqData limsShallowSeqData = shallowSeqPuritySample.get(barcode);

        assertEquals(barcode, limsShallowSeqData.sampleBarcode());
        assertEquals("SAMP01011234T", limsShallowSeqData.sampleId());
        assertEquals("0.19", limsShallowSeqData.purityShallowSeq());
        assertTrue(limsShallowSeqData.hasReliablePurity());
        assertTrue(limsShallowSeqData.hasReliableQuality());
    }
}