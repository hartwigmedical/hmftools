package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.sampleBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class BiopsyMatcherTest {
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private final static LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private final static LocalDate AUG2015 = LocalDate.parse("2015-08-01");
    private final static LocalDate SEP2015 = LocalDate.parse("2015-09-01");
    private final static LocalDate NOV2015 = LocalDate.parse("2015-11-01");

    private final static SampleData LIMS_SAMPLE_JUL = sampleBuilder(JUL2015).build();
    private final static SampleData LIMS_SAMPLE_AUG = sampleBuilder(AUG2015).build();
    private final static SampleData LIMS_SAMPLE_SEP = sampleBuilder(SEP2015).build();
    private final static SampleData LIMS_SAMPLE_NOV = sampleBuilder(NOV2015).build();

    private final static SampleData LIMS_ARRIVED_NOV_SAMPLED_MAR = sampleBuilder(NOV2015).samplingDate(MAR2015).build();

    private final static BiopsyData ECRF_BIOPSY_JAN = biopsyBuilder().date(JAN2015).build();
    private final static BiopsyData ECRF_BIOPSY_FEB = biopsyBuilder().date(FEB2015).build();
    private final static BiopsyData ECRF_BIOPSY_MAR = biopsyBuilder().date(MAR2015).build();
    private final static BiopsyData ECRF_BIOPSY_MAR_NOT_EVALUABLE = biopsyBuilder().date(MAR2015).biopsyEvaluable("no").build();
    private final static BiopsyData ECRF_BIOPSY_JUL = biopsyBuilder().date(JUL2015).build();
    private final static BiopsyData ECRF_BIOPSY_SEP = biopsyBuilder().date(SEP2015).build();

    private final static BiopsyData BIOPSY_NULL = biopsyBuilder().build();

    // MIVO:    ---biopsy(jul)/sample(jul)---
    @Test
    public void biopsyAndSampleSameDateYieldsMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_JUL);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, LIMS_SAMPLE_JUL.sampleId());
    }

    // MIVO:    ---biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void biopsyBeforeSampleWithinThresholdYieldsMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, LIMS_SAMPLE_JUL.sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-x-x-x-x-sample(jul)---
    @Test
    public void biopsyBeforeSampleOutsideThresholdYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_JAN);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null);
    }

    // MIVO:    ---sample(aug)-biopsy(sep)---
    @Test
    public void biopsyTakenAfterSampleArrivedYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_AUG);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null);
    }

    // MIVO:    ---biopsy(feb)-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesBeforeSampleWithinThresholdYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_FEB, ECRF_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    // MIVO:    ---biopsy(feb)-x-x-null(?)-x-sample(jul)---
    @Test
    public void twoBiopsiesSecondNullBeforeSampleWithinThresholdYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_FEB, BIOPSY_NULL);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    // MIVO:    ---null(?)-biopsy(feb)-x-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesFirstNullBeforeSampleWithinThresholdYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(BIOPSY_NULL, ECRF_BIOPSY_FEB);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    // MIVO:    ---biopsy(jan)-x-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesBeforeSampleOneWithinThresholdYieldsMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_JAN, ECRF_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, LIMS_SAMPLE_JUL.sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesBeforeSample1WithinThresholdOutOfOrderYieldsMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_MAR, ECRF_BIOPSY_JAN);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, LIMS_SAMPLE_JUL.sampleId());
    }

    // MIVO:    ---biopsy(jan)-biopsy(feb)-x-x-x-x-x-x-sample(sep)---
    @Test
    public void twoBiopsiesBeforeSampleOutsideThresholdYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_SEP);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_JAN, ECRF_BIOPSY_FEB);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    // MIVO:    ---biopsy(feb)-x-x-x-x-sample(jul)-sample(aug)---
    @Test
    public void biopsyBefore2SamplesWithinThresholdYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL, LIMS_SAMPLE_AUG);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_FEB);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null);
    }

    // MIVO:    ---biopsy(mar)-x-x-x-sample(jul)-x-biopsy(sep)---
    @Test
    public void oneSampleBetween2BiopsiesWithinThresholdYieldsMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_MAR, ECRF_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, LIMS_SAMPLE_JUL.sampleId(), null);
    }

    // MIVO:    ---biopsy(feb)-biopsy(mar)-x-x-x-x-sample(aug)-sample(sep)---
    @Test
    public void twoSamplesAfterTwoBiopsiesYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_AUG, LIMS_SAMPLE_SEP);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_FEB, ECRF_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    // MIVO:    ---biopsy(mar)-x-x-x-x-sample(aug)-biopsy(sep)-x-sample(nov)---
    @Test
    public void twoSamplesBetweenTwoBiopsiesWithinThresholdYieldsMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_AUG, LIMS_SAMPLE_NOV);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_MAR, ECRF_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, LIMS_SAMPLE_AUG.sampleId(), LIMS_SAMPLE_NOV.sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-x-x-x-x-x-sample(aug)-biopsy(sep)-x-sample(nov)---
    @Test
    public void twoSamplesBetweenTwoBiopsiesOneOutsideThresholdYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_AUG, LIMS_SAMPLE_NOV);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_JAN, ECRF_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    // KODU: --biopsy(mar)-x-x-x-x-sample(nov)--
    @Test
    public void biopsyArrivedInNovButSamplingDateKnownYieldsMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_ARRIVED_NOV_SAMPLED_MAR);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, LIMS_ARRIVED_NOV_SAMPLED_MAR.sampleId());
    }

    // KODU: --biopsy(jan)-x-sample(mar)-x-x-sample_arrived(nov)--
    @Test
    public void biopsyFromJanArrivedInNovButSamplingDateInMarYieldsNoMatch() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_ARRIVED_NOV_SAMPLED_MAR);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_JAN);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null);
    }

    @Test
    public void noMatchWhenBiopsyIsNotEvaluable() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyData> clinicalBiopsies = Lists.newArrayList(ECRF_BIOPSY_MAR_NOT_EVALUABLE);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null);
    }

    private static void runMatcherAndVerify(@NotNull List<SampleData> sequencedBiopsies, @NotNull List<BiopsyData> clinicalBiopsies,
            @Nullable String... expectedSampleIds) {
        final List<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples("patient", sequencedBiopsies, clinicalBiopsies).values();
        Collections.sort(matchedBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        if (expectedSampleIds == null) {
            assertEquals(null, matchedBiopsies.get(0).sampleId());
        } else {
            for (int i = 0; i < expectedSampleIds.length; i++) {
                assertTrue(matchedBiopsies.size() >= i);
                assertEquals(expectedSampleIds[i], matchedBiopsies.get(i).sampleId());
            }
        }
    }
}
