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
import org.junit.Ignore;
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
        final String label = "WIDE";
        final String projectName = "projectX";
        final String contactEmail = "henk@hmf.nl";
        final String contactName = "henk";

        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(SAMPLE)
                .samplingDateString(samplingDate)
                .arrivalDateString(arrivalDate)
                .dnaConcentration(dnaConcentration)
                .tumorPercentageString(tumorPercentage)
                .primaryTumor(primaryTumor)
                .labSopVersions(labSopVersions)
                .labRemarks(labRemarks)
                .labelSample(label)
                .projectName(projectName)
                .submission(SUBMISSION)
                .build();

        final LimsJsonSubmissionData submissionData = ImmutableLimsJsonSubmissionData.builder()
                .submission(SUBMISSION)
                .contactEmail(contactEmail)
                .contactName(contactName)
                .build();

        final Lims lims = buildTestLimsWithSampleAndSubmission(sampleData, submissionData);

        assertEquals(1, lims.sampleCount());

        assertEquals(contactEmail, lims.contactEmail(SAMPLE));
        assertEquals(contactName, lims.contactName(SAMPLE));
        assertNull(lims.patientNumber(SAMPLE));
        assertEquals(label, lims.labelSample(SAMPLE));
        assertEquals(projectName, lims.projectNameDVO(SAMPLE));
        assertEquals(LimsTestUtil.toDate(arrivalDate), lims.arrivalDateForSample(SAMPLE));
        assertEquals(LimsTestUtil.toDate(samplingDate), lims.samplingDateForSample(SAMPLE));

        Integer dnaAmount = lims.dnaNanogramsForSample(SAMPLE);
        assertNotNull(dnaAmount);
        assertEquals(500L, (int) dnaAmount);

        String tumorPerc = lims.tumorPercentageForSample(SAMPLE);
        assertEquals("40%", tumorPerc);

        assertEquals(primaryTumor, lims.primaryTumorForSample(SAMPLE));
        assertEquals(labSopVersions, lims.labProceduresForSample(SAMPLE));
    }

    @Test
    public void worksForNonExistingSamplesAndSubmissions() {
        Lims lims = LimsFactory.empty();

        assertEquals("N/A", lims.contactName("DoesNotExist"));
        assertEquals("N/A", lims.contactEmail("DoesNotExist"));
        assertNull(lims.patientNumber("DoesNotExist"));
        assertEquals("N/A", lims.labelSample("DoesNotExist"));
        assertEquals("N/A", lims.projectNameDVO("DoesNotExist"));
        assertNull(lims.arrivalDateForSample("DoesNotExist"));
        assertNull(lims.samplingDateForSample("DoesNotExist"));
        assertNull(lims.dnaNanogramsForSample("DoesNotExist"));
        assertEquals("N/A", lims.tumorPercentageForSample("DoesNotExist"));
        assertEquals("N/A", lims.primaryTumorForSample("DoesNotExist"));
        assertEquals("N/A", lims.labProceduresForSample("DoesNotExist"));
    }

    @Test
    public void fallBackOnPreLIMSArrivalDateWorks() {
        final LocalDate date = LimsTestUtil.toDate("2017-10-03");

        final Lims lims = buildTestLimsWithPreLIMSArrivalDateForSample(SAMPLE, date);

        assertEquals(date, lims.arrivalDateForSample(SAMPLE));
    }

    @Test
    public void invalidDataYieldsNullOrNA() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder()
                .sampleId(SAMPLE)
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

        assertEquals("N/A", lims.contactEmail(SAMPLE));
        assertEquals("N/A", lims.contactName(SAMPLE));
        assertNull(lims.arrivalDateForSample(SAMPLE));
        assertNull(lims.samplingDateForSample(SAMPLE));
        assertNull(lims.dnaNanogramsForSample(SAMPLE));
        assertEquals("N/A", lims.tumorPercentageForSample(SAMPLE));
        assertEquals("N/A", lims.labProceduresForSample(SAMPLE));
    }

    @Test
    public void noTumorPercentageForCORE() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder()
                .sampleId(SAMPLE)
                .labelSample("CORE")
                .build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("not determined", lims.tumorPercentageForSample(SAMPLE));
    }

    @Test
    public void noTumorPercentageForShallowSeq() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder()
                .sampleId(SAMPLE)
                .labRemarks("ShallowSeq")
                .build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("not determined", lims.tumorPercentageForSample(SAMPLE));
    }

    //TODO: tests for shallow seq purity
    @Ignore
    public void ShallowSeqPurityError() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder()
                .sampleId(SAMPLE)
                .labRemarks("ShallowSeq")
                .build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("N/A", lims.purityShallowSeq(SAMPLE));
    }

    @Ignore
    public void ShallowSeqPurity() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder()
                .sampleId("CPCT02990001T")
                .labRemarks("ShallowSeq")
                .build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("N/A", lims.purityShallowSeq("CPCT02990001T"));
    }

    @Ignore
    public void ShallowSeqPathology() {
        final LimsJsonSampleData sampleData = createLimsSampleDataBuilder()
                .sampleId(SAMPLE)
                .labRemarks("ShallowSeq")
                .build();

        Lims lims = buildTestLimsWithSample(sampleData);
        assertEquals("N/A", lims.purityShallowSeq(SAMPLE));
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