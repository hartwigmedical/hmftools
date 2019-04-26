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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LimsTest {

    private static final String SAMPLE = "CPCT02991111T";
    private static final String SUBMISSION = "ABCDEF123";

    @Test
    public void canReadProperlyDefinedSample() {
        final String patientId = "CPCT02991111";
        final String arrivalDate = "2017-05-01";
        final String samplingDate = "2017-04-15";
        final String dnaConcentration = "10";
        final String purityShallowSeq = "0.71";
        final String primaryTumor = "Prostate";
        final String labSopVersions = "PREP1V2-QC1V2-SEQ1V2";
        final String projectName = "projectX";
        final String requesterEmail = "henk@hmf.nl";
        final String requesterName = "henk";
        final String refBarcode = "A123";
        final String tumorBarcode = "T456";
        final String labRemarks = "CPCT WIDE project";
        final String hospitalPatientId = "Henkie";
        final String hospitalPathologySampleId = "Henkie's sample";

        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE)
                .patientId(patientId)
                .arrivalDate(arrivalDate)
                .samplingDate(samplingDate)
                .dnaConcentration(dnaConcentration)
                .primaryTumor(primaryTumor)
                .labSopVersions(labSopVersions)
                .labRemarks(labRemarks)
                .submission(SUBMISSION)
                .refBarcode(refBarcode)
                .tumorBarcode(tumorBarcode)
                .requesterEmail(requesterEmail)
                .requesterName(requesterName)
                .shallowSeq("1")
                .germlineReportingChoice(Strings.EMPTY)
                .hospitalPatientId(hospitalPatientId)
                .hospitalPathologySampleId(hospitalPathologySampleId)
                .build();

        final LimsJsonSubmissionData submissionData =
                ImmutableLimsJsonSubmissionData.builder().submission(SUBMISSION).projectName(projectName).build();

        final LimsShallowSeqData shallowSeqData =
                ImmutableLimsShallowSeqData.builder().sampleId(SAMPLE).purityShallowSeq(purityShallowSeq).build();

        final Lims lims = buildFullTestLims(sampleData, submissionData, shallowSeqData);

        assertEquals(1, lims.sampleCount());

        assertEquals(patientId, lims.patientId(SAMPLE));
        assertEquals(refBarcode, lims.refBarcode(SAMPLE));
        assertEquals(tumorBarcode, lims.tumorBarcode(SAMPLE));
        assertEquals(LimsTestUtil.toDate(samplingDate), lims.samplingDate(SAMPLE));
        assertEquals(LimsTestUtil.toDate(arrivalDate), lims.arrivalDate(SAMPLE));

        assertEquals(SUBMISSION, lims.submissionId(SAMPLE));
        assertEquals(projectName, lims.projectName(SAMPLE));
        assertEquals(requesterEmail, lims.requesterEmail(SAMPLE));
        assertEquals(requesterName, lims.requesterName(SAMPLE));

        Integer dnaAmount = lims.dnaNanograms(SAMPLE);
        assertNotNull(dnaAmount);
        assertEquals(500L, (int) dnaAmount);

        assertEquals("71%", lims.purityShallowSeq(SAMPLE));
        assertEquals(Lims.NOT_DETERMINED_STRING, lims.pathologyTumorPercentage(SAMPLE));
        assertEquals(primaryTumor, lims.primaryTumor(SAMPLE));
        assertEquals(labSopVersions, lims.labProcedures(SAMPLE));

        assertEquals(hospitalPatientId, lims.hospitalPatientId(SAMPLE));
        assertEquals(hospitalPathologySampleId, lims.hospitalPathologySampleId(SAMPLE));
        assertEquals(LimsGermlineReportingChoice.UNKNOWN, lims.germlineReportingChoice(SAMPLE));
    }

    @Test
    public void worksForNonExistingSamplesAndSubmissions() {
        Lims lims = LimsFactory.empty();
        String doesNotExistSample = "DoesNotExist";

        assertNull(lims.arrivalDate(doesNotExistSample));
        assertNull(lims.samplingDate(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.refBarcode(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.tumorBarcode(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.requesterName(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.requesterEmail(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.projectName(doesNotExistSample));
        assertNull(lims.dnaNanograms(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.pathologyTumorPercentage(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.purityShallowSeq(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.primaryTumor(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.labProcedures(doesNotExistSample));
        assertEquals(Lims.NOT_KNOWN_STRING, lims.hospitalPatientId(doesNotExistSample));
        assertEquals(Lims.NOT_KNOWN_STRING, lims.hospitalPathologySampleId(doesNotExistSample));
        assertEquals(LimsGermlineReportingChoice.UNKNOWN, lims.germlineReportingChoice(doesNotExistSample));
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
                .arrivalDate("IsNotADate")
                .samplingDate(null)
                .dnaConcentration("IsNotADNAConcentration")
                .pathologyTumorPercentage("IsNotANumber")
                .labSopVersions("anything")
                .build();

        final Lims lims = buildTestLimsWithSample(sampleData);

        assertEquals(1, lims.sampleCount());

        assertNull(lims.arrivalDate(SAMPLE));
        assertNull(lims.samplingDate(SAMPLE));
        assertNull(lims.dnaNanograms(SAMPLE));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.pathologyTumorPercentage(SAMPLE));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.labProcedures(SAMPLE));
    }

    @Test
    public void noPathologyTumorPercentageDeterminedForShallowSeq() {
        final LimsJsonSampleData sampleData1 = createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq("1").build();

        Lims lims1 = buildTestLimsWithSample(sampleData1);
        assertEquals(Lims.NOT_DETERMINED_STRING, lims1.pathologyTumorPercentage(SAMPLE));

        final LimsJsonSampleData sampleData2 = createLimsSampleDataBuilder().sampleId(SAMPLE).labRemarks("ShallowSeq").build();

        Lims lims2 = buildTestLimsWithSample(sampleData2);
        assertEquals(Lims.NOT_DETERMINED_STRING, lims2.pathologyTumorPercentage(SAMPLE));
    }

    @Test
    public void missingOrMalformedShallowSeqDataForSampleYieldsNA() {
        final LimsJsonSampleData sampleData1 = createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq("1").build();

        Lims lims1 = buildTestLimsWithSample(sampleData1);
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims1.purityShallowSeq(SAMPLE));
        assertEquals(Lims.NOT_DETERMINED_STRING, lims1.pathologyTumorPercentage(SAMPLE));

        final LimsJsonSampleData sampleData2 = createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq("1").build();
        Lims lims2 = buildTestLimsWithSampleAndShallowSeq(sampleData2, "NotANumber");
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims2.purityShallowSeq(SAMPLE));
        assertEquals(Lims.NOT_DETERMINED_STRING, lims2.pathologyTumorPercentage(SAMPLE));
    }

    @Test
    public void canRetrievePathologyPercentageForSample() {
        final LimsJsonSampleData sampleData1 =
                createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq("0").pathologyTumorPercentage("70").build();

        Lims lims1 = buildTestLimsWithSample(sampleData1);
        assertEquals(Lims.NOT_DETERMINED_STRING, lims1.purityShallowSeq(SAMPLE));
        assertEquals("70%", lims1.pathologyTumorPercentage(SAMPLE));

        final LimsJsonSampleData sampleData2 =
                createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq("0").pathologyTumorPercentage("NotANumber").build();

        Lims lims2 = buildTestLimsWithSample(sampleData2);
        assertEquals(Lims.NOT_DETERMINED_STRING, lims2.purityShallowSeq(SAMPLE));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims2.pathologyTumorPercentage(SAMPLE));
    }

    @Test
    public void canRetrieveShallowSeqPurityForSample() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq("1").build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "0.2");
        assertEquals("20%", lims.purityShallowSeq(SAMPLE));
        assertEquals(Lims.NOT_DETERMINED_STRING, lims.pathologyTumorPercentage(SAMPLE));
    }

    @Test
    public void canRetrieveShallowSeqBelowDetectionLimitForSample() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE).shallowSeq("1").build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "below detection threshold");
        assertEquals("below detection threshold", lims.purityShallowSeq(SAMPLE));
        assertEquals(Lims.NOT_DETERMINED_STRING, lims.pathologyTumorPercentage(SAMPLE));
    }

    @NotNull
    private static Lims buildFullTestLims(@NotNull final LimsJsonSampleData sampleData,
            @NotNull final LimsJsonSubmissionData submissionData, @NotNull LimsShallowSeqData shallowSeqData) {
        Map<String, LimsJsonSampleData> dataPerSample = Maps.newHashMap();
        dataPerSample.put(sampleData.sampleId(), sampleData);

        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        dataPerSubmission.put(submissionData.submission(), submissionData);

        Map<String, LocalDate> preLIMSArrivalDates = Maps.newHashMap();
        Set<String> samplesWithSamplingDates = Sets.newHashSet();

        Map<String, LimsShallowSeqData> shallowSeqDataPerSample = Maps.newHashMap();
        shallowSeqDataPerSample.put(shallowSeqData.sampleId(), shallowSeqData);

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
    private static Lims buildTestLimsWithSampleAndShallowSeq(@NotNull final LimsJsonSampleData sampleData,
            @NotNull String shallowSeqPurity) {
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