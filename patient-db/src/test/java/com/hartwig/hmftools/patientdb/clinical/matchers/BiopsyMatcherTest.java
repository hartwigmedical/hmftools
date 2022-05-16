package com.hartwig.hmftools.patientdb.clinical.matchers;

import static com.hartwig.hmftools.patientdb.clinical.datamodel.DatamodelTestFactory.biopsyBuilder;
import static com.hartwig.hmftools.patientdb.clinical.datamodel.DatamodelTestFactory.sampleBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class BiopsyMatcherTest {

    private static final LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private static final LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private static final LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private static final LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private static final LocalDate AUG2015 = LocalDate.parse("2015-08-01");
    private static final LocalDate SEP2015 = LocalDate.parse("2015-09-01");
    private static final LocalDate NOV2015 = LocalDate.parse("2015-11-01");

    private static final SampleData SEQUENCED_BIOPSY_JUL = sampleBuilder(JUL2015).build();
    private static final SampleData SEQUENCED_BIOPSY_AUG = sampleBuilder(AUG2015).build();
    private static final SampleData SEQUENCED_BIOPSY_SEP = sampleBuilder(SEP2015).build();
    private static final SampleData SEQUENCED_BIOPSY_NOV = sampleBuilder(NOV2015).build();

    private static final SampleData SEQUENCED_BIOPSY_ARRIVED_NOV_SAMPLED_MAR = sampleBuilder(NOV2015).samplingDate(MAR2015).build();

    private static final BiopsyData CLINICAL_BIOPSY_JAN = biopsyBuilder().date(JAN2015).build();
    private static final BiopsyData CLINICAL_BIOPSY_FEB = biopsyBuilder().date(FEB2015).build();
    private static final BiopsyData CLINICAL_BIOPSY_MAR = biopsyBuilder().date(MAR2015).build();
    private static final BiopsyData CLINICAL_BIOPSY_MAR_NOT_EVALUABLE = biopsyBuilder().date(MAR2015).biopsyEvaluable("no").build();
    private static final BiopsyData CLINICAL_BIOPSY_JUL = biopsyBuilder().date(JUL2015).build();
    private static final BiopsyData CLINICAL_BIOPSY_SEP = biopsyBuilder().date(SEP2015).build();

    private static final BiopsyData CLINICAL_BIOPSY_NULL = biopsyBuilder().build();

    //    ---biopsy(jul)/sample(jul)---
    @Test
    public void biopsyAndSampleSameDateYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_JUL);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, SEQUENCED_BIOPSY_JUL.sampleId());
    }

    //    ---biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void biopsyBeforeSampleWithinThresholdYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, SEQUENCED_BIOPSY_JUL.sampleId());
    }

    //    ---biopsy(jan)-x-x-x-x-x-sample(jul)---
    @Test
    public void biopsyBeforeSampleOutsideThresholdYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_JAN);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, (String) null);
    }

    //    ---sample(aug)-biopsy(sep)---
    @Test
    public void biopsyTakenAfterSampleArrivedYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_AUG);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, (String) null);
    }

    //    ---biopsy(feb)-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesBeforeSampleWithinThresholdYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_FEB, CLINICAL_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    //    ---biopsy(feb)-x-x-null(?)-x-sample(jul)---
    @Test
    public void twoBiopsiesSecondNullBeforeSampleWithinThresholdYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_FEB, CLINICAL_BIOPSY_NULL);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    //   ---null(?)-biopsy(feb)-x-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesFirstNullBeforeSampleWithinThresholdYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_NULL, CLINICAL_BIOPSY_FEB);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    //    ---biopsy(jan)-x-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesBeforeSampleOneWithinThresholdYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_JAN, CLINICAL_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, SEQUENCED_BIOPSY_JUL.sampleId());
    }

    //    ---biopsy(jan)-x-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void twoBiopsiesBeforeSample1WithinThresholdOutOfOrderYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_MAR, CLINICAL_BIOPSY_JAN);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, SEQUENCED_BIOPSY_JUL.sampleId());
    }

    //    ---biopsy(jan)-biopsy(feb)-x-x-x-x-x-x-sample(sep)---
    @Test
    public void twoBiopsiesBeforeSampleOutsideThresholdYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_SEP);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_JAN, CLINICAL_BIOPSY_FEB);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    //    ---biopsy(feb)-x-x-x-x-sample(jul)-sample(aug)---
    @Test
    public void biopsyBefore2SamplesWithinThresholdYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL, SEQUENCED_BIOPSY_AUG);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_FEB);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, SEQUENCED_BIOPSY_JUL.sampleId());
    }

    //    ---biopsy(mar)-x-x-x-sample(jul)-x-biopsy(sep)---
    @Test
    public void oneSampleBetween2BiopsiesWithinThresholdYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_MAR, CLINICAL_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, SEQUENCED_BIOPSY_JUL.sampleId(), null);
    }

    //    ---biopsy(feb)-biopsy(mar)-x-x-x-x-sample(aug)-sample(sep)---
    @Test
    public void twoSamplesAfterTwoBiopsiesYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_SEP, SEQUENCED_BIOPSY_NOV);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_FEB, CLINICAL_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, null);
    }

    //    ---biopsy(mar)-x-x-x-x-sample(aug)-biopsy(sep)-x-sample(nov)---
    @Test
    public void twoSamplesBetweenTwoBiopsiesWithinThresholdYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_AUG, SEQUENCED_BIOPSY_NOV);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_MAR, CLINICAL_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, SEQUENCED_BIOPSY_AUG.sampleId(), SEQUENCED_BIOPSY_NOV.sampleId());
    }

    //    ---biopsy(jan)-x-x-x-x-x-x-sample(aug)-biopsy(sep)-x-sample(nov)---
    @Test
    public void twoSamplesBetweenTwoBiopsiesOneOutsideThresholdYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_AUG, SEQUENCED_BIOPSY_NOV);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_JAN, CLINICAL_BIOPSY_SEP);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, null, SEQUENCED_BIOPSY_NOV.sampleId());
    }

    // --biopsy(mar)-x-x-x-x-sample(nov)--
    @Test
    public void biopsyArrivedInNovButSamplingDateKnownYieldsMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_ARRIVED_NOV_SAMPLED_MAR);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_MAR);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, SEQUENCED_BIOPSY_ARRIVED_NOV_SAMPLED_MAR.sampleId());
    }

    // --biopsy(jan)-x-sample(mar)-x-x-sample_arrived(nov)--
    @Test
    public void biopsyFromJanArrivedInNovButSamplingDateInMarYieldsNoMatch() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_ARRIVED_NOV_SAMPLED_MAR);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_JAN);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, (String) null);
    }

    @Test
    public void noMatchWhenBiopsyIsNotEvaluable() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList(CLINICAL_BIOPSY_MAR_NOT_EVALUABLE);

        runMatcherAndVerify(sequencedBiopsies, clinicalBiopsies, (String) null);
    }

    @Test
    public void createFindingWhenNotEnoughClinicalBiopsies() {
        List<SampleData> sequencedBiopsies = Lists.newArrayList(SEQUENCED_BIOPSY_JUL);
        List<BiopsyData> clinicalBiopsies = Lists.newArrayList();

        MatchResult<BiopsyData> matchedBiopsies = BiopsyMatcher.matchBiopsiesToTumorSamples("patient", sequencedBiopsies, clinicalBiopsies);

        assertEquals(2, matchedBiopsies.findings().size());
        assertTrue(matchedBiopsies.findings().get(0).message().toLowerCase().contains("not enough clinical biopsy forms to match"));
    }

    private static void runMatcherAndVerify(@NotNull List<SampleData> sequencedBiopsies, @NotNull List<BiopsyData> clinicalBiopsies,
            @Nullable String... expectedSampleIds) {
        List<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples("patient", sequencedBiopsies, clinicalBiopsies).values();
        assertEquals(clinicalBiopsies.size(), matchedBiopsies.size());

        if (expectedSampleIds == null) {
            assertNull(matchedBiopsies.get(0).sampleId());
        } else {
            Collections.sort(matchedBiopsies);
            for (int i = 0; i < expectedSampleIds.length; i++) {
                assertTrue(matchedBiopsies.size() >= i);
                assertEquals(expectedSampleIds[i], matchedBiopsies.get(i).sampleId());
            }
        }
    }
}
