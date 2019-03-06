package com.hartwig.hmftools.common.lims;

import static com.hartwig.hmftools.common.lims.LimsTestUtil.createLimsSampleDataBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LimsTest {

    private static final String SAMPLE = "CPCT02991111T";
    private static final String SUBMISSION = "ABCDEF123";

    @Test
    public void canReadProperlyDefinedSample() {
        final String arrivalDate = "2017-05-01";
        final String samplingDate = "2017-04-15";
        final String dnaConcentration = "10";
        final String tumorPercentage = "40";
        final String primaryTumor = "Prostate";
        final String labSopVersions = "PREP1V2-QC1V2-SEQ1V2";
        final String labRemarks = "CPCT WIDE project";
        final String projectName = "projectX";
        final String contactEmail = "henk@hmf.nl";
        final String contactName = "henk";
        final String refBarcode = "A123";
        final String tumorBarcode = "not determinded";
        final String patientId = "CPCT02991111";


        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE)
                .samplingDateString(samplingDate)
                .arrivalDateString(arrivalDate)
                .dnaConcentration(dnaConcentration)
                .tumorPercentageString(tumorPercentage)
                .primaryTumor(primaryTumor)
                .labSopVersions(labSopVersions)
                .labRemarks(labRemarks)
                .projectName(projectName)
                .submissionSamples(SUBMISSION)
                .refBarcodeId(refBarcode)
                .tumorBarcodeId(tumorBarcode)
                .patientId(patientId)
                .requesterEmail(contactEmail)
                .requesterName(contactName)
                .shallowSeq(0)
                .build();

        final LimsJsonSubmissionData submissionData = ImmutableLimsJsonSubmissionData.builder()
                .submission(SUBMISSION)
                .projectName("projectX")
                .build();

        final Lims lims = buildTestLimsWithSampleAndSubmission(sampleData, submissionData);

        assertEquals(1, lims.sampleCount());

        assertEquals(contactEmail, lims.contactEmails(SAMPLE));
        assertEquals(contactName, lims.contactNames(SAMPLE));
        assertNull(lims.hospitalPatientId(SAMPLE));
        assertEquals(projectName, lims.projectName(SAMPLE));
        assertEquals(LimsTestUtil.toDate(arrivalDate), lims.arrivalDate(SAMPLE));
        assertEquals(LimsTestUtil.toDate(samplingDate), lims.samplingDate(SAMPLE));
        assertEquals(refBarcode, lims.barcodeReference(SAMPLE));
        assertEquals(tumorBarcode, lims.barcodeTumor(SAMPLE));

        Integer dnaAmount = lims.dnaNanograms(SAMPLE);
        assertNotNull(dnaAmount);
        assertEquals(500L, (int) dnaAmount);

        assertEquals("40%", lims.pathologyTumorPercentage(SAMPLE));
        assertEquals(primaryTumor, lims.primaryTumor(SAMPLE));
        assertEquals(labSopVersions, lims.labProcedures(SAMPLE));
    }

    @Test
    public void worksForNonExistingSamplesAndSubmissions() {
        Lims lims = LimsFactory.empty();

        assertEquals("N/A", lims.contactNames("DoesNotExist"));
        assertEquals("N/A", lims.contactEmails("DoesNotExist"));
        assertNull(lims.hospitalPatientId("DoesNotExist"));
        assertEquals("N/A", lims.projectName("DoesNotExist"));
        assertNull(lims.arrivalDate("DoesNotExist"));
        assertNull(lims.samplingDate("DoesNotExist"));
        assertNull(lims.dnaNanograms("DoesNotExist"));
        assertEquals("N/A", lims.pathologyTumorPercentage("DoesNotExist"));
        assertEquals("N/A", lims.purityShallowSeq("DoesNotExist"));
        assertEquals("N/A", lims.primaryTumor("DoesNotExist"));
        assertEquals("N/A", lims.labProcedures("DoesNotExist"));
    }

    @Test
    public void fallBackOnPreLIMSArrivalDateWorks() {
        final LocalDate date = LimsTestUtil.toDate("2017-10-03");

        final Lims lims = buildTestLimsWithPreLIMSArrivalDateForSample(SAMPLE, date);

        assertEquals(date, lims.arrivalDate(SAMPLE));
    }

    @Test
    public void invalidDataYieldsNullOrNA() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE)
                .samplingDateString("IsNotADate")
                .arrivalDateString("IsNotADate")
                .dnaConcentration("IsNotADNAConcentration")
                .tumorPercentageString("IsNotANumber")
                .primaryTumor("anything")
                .labSopVersions("anything")
                .labRemarks("anything")
                .build();

        final Lims lims = buildTestLimsWithSample(sampleData);

        assertEquals(1, lims.sampleCount());

        assertEquals("", lims.contactEmails(SAMPLE));
        assertEquals("", lims.contactNames(SAMPLE));
        assertNull(lims.arrivalDate(SAMPLE));
        assertNull(lims.samplingDate(SAMPLE));
        assertNull(lims.dnaNanograms(SAMPLE));
        assertEquals("N/A", lims.pathologyTumorPercentage(SAMPLE));
        assertEquals("not determined", lims.purityShallowSeq(SAMPLE));
        assertEquals("N/A", lims.labProcedures(SAMPLE));
    }

    @Test
    public void noPathologyTumorPercentageDeterminedForCORE() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq(1).build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("not determined", lims.pathologyTumorPercentage(SAMPLE));
    }

    @Test
    public void noPathologyTumorPercentageForShallowSeq() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE).labRemarks("ShallowSeq").build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("not determined", lims.pathologyTumorPercentage(SAMPLE));
    }

    @Test
    public void shallowSeqPurityErrorForSample() {
        final LimsJsonSampleData sampleData =
                createLimsSampleDataBuilder().sampleId(SAMPLE).labelSample("").labRemarks("ShallowSeq").build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("N/A", lims.purityShallowSeq(SAMPLE));
    }

    @Test
    public void canRetrievePathologyPercentageForSample() {
        final LimsJsonSampleData sampleData =
                createLimsSampleDataBuilder().sampleId(SAMPLE).labRemarks("").shallowSeq(0).tumorPercentageString("70").build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("not determined", lims.purityShallowSeq(SAMPLE));
        assertEquals("70%", lims.pathologyTumorPercentage(SAMPLE));
    }

    @Test
    public void canRetrieveShallowSeqPurityForCPCTSample() {
        final LimsJsonSampleData sampleData =
                createLimsSampleDataBuilder().sampleId("CPCT02990001T").labRemarks("ShallowSeq").build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "0.2");
        assertEquals("20%", lims.purityShallowSeq("CPCT02990001T"));
    }

    @Test
    public void canRetrieveShallowSeqPurityForCORESample() {
        final LimsJsonSampleData sampleData =
                createLimsSampleDataBuilder().sampleId("CORE00000001T").labRemarks("ShallowSeq").labelSample("CORE").build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "0.3");
        assertEquals("30%", lims.purityShallowSeq("CORE00000001T"));
    }

    @Test
    public void canRetrieveShallowSeqBelowDetectionLimitForCPCTSample() {
        final LimsJsonSampleData sampleData =
                createLimsSampleDataBuilder().sampleId("CPCT02990003T").labelSample("").labRemarks("ShallowSeq").build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "below detection threshold");
        assertEquals("below detection threshold", lims.purityShallowSeq("CPCT02990003T"));
    }

    @Test
    public void canRetrieveShallowSeqBelowDetectionLimitForCORESample() {
        final LimsJsonSampleData sampleData =
                createLimsSampleDataBuilder().sampleId("CPCT02990003T").labelSample("").labRemarks("ShallowSeq").build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "below detection threshold");
        assertEquals("below detection threshold", lims.purityShallowSeq("CPCT02990003T"));
    }

    @NotNull
    private static Lims buildTestLimsWithSampleAndSubmission(@NotNull final LimsJsonSampleData sampleData,
            @NotNull final LimsJsonSubmissionData submissionData) {
        Map<String, LimsJsonSampleData> dataPerSample = Maps.newHashMap();
        dataPerSample.put(sampleData.sampleId(), sampleData);
        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        dataPerSubmission.put(submissionData.submission(), submissionData);
        Map<String, LocalDate> preLIMSArrivalDates = Maps.newHashMap();
        Set<String> samplesWithSamplingDates = Sets.newHashSet();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSample = Maps.newHashMap();

        return new Lims(dataPerSample, dataPerSubmission, preLIMSArrivalDates, samplesWithSamplingDates, shallowSeqDataPerSample);
    }

    @NotNull
    private static Lims buildTestLimsWithSample(@NotNull final LimsJsonSampleData sampleData) {
        Map<String, LimsJsonSampleData> dataPerSample = Maps.newHashMap();
        dataPerSample.put(sampleData.sampleId(), sampleData);
        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        Map<String, LocalDate> preLIMSArrivalDates = Maps.newHashMap();
        Set<String> samplesWithSamplingDates = Sets.newHashSet();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSample = Maps.newHashMap();

        return new Lims(dataPerSample, dataPerSubmission, preLIMSArrivalDates, samplesWithSamplingDates, shallowSeqDataPerSample);
    }

    @NotNull
    private static Lims buildTestLimsWithSampleAndShallowSeq(@NotNull final LimsJsonSampleData sampleData, String shallowSeqPurity) {
        Map<String, LimsJsonSampleData> dataPerSample = Maps.newHashMap();
        dataPerSample.put(sampleData.sampleId(), sampleData);
        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        Map<String, LocalDate> preLIMSArrivalDates = Maps.newHashMap();
        Set<String> samplesWithSamplingDates = Sets.newHashSet();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSample = Maps.newHashMap();

        shallowSeqDataPerSample.put(sampleData.sampleId(), ImmutableLimsShallowSeqData.of(sampleData.sampleId(), shallowSeqPurity));

        return new Lims(dataPerSample, dataPerSubmission, preLIMSArrivalDates, samplesWithSamplingDates, shallowSeqDataPerSample);
    }

    @NotNull
    private static Lims buildTestLimsWithPreLIMSArrivalDateForSample(@NotNull final String sample, @NotNull final LocalDate date) {
        final Map<String, LimsJsonSampleData> dataPerSample = Maps.newHashMap();
        final Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        final Map<String, LocalDate> preLIMSArrivalDates = Maps.newHashMap();
        preLIMSArrivalDates.put(sample, date);

        Set<String> samplesWithSamplingDates = Sets.newHashSet();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSample = Maps.newHashMap();

        return new Lims(dataPerSample, dataPerSubmission, preLIMSArrivalDates, samplesWithSamplingDates, shallowSeqDataPerSample);
    }
}