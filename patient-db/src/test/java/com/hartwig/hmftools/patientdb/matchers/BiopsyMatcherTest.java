package com.hartwig.hmftools.patientdb.matchers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyData;
import com.hartwig.hmftools.patientdb.data.ImmutableSampleData;
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

    private final static SampleData LIMS_SAMPLE_JUL = ImmutableSampleData.of("jul-sample", JUL2015, null, 0D);
    private final static SampleData LIMS_SAMPLE_AUG = ImmutableSampleData.of("aug-sample", AUG2015, null, 0D);
    private final static SampleData LIMS_SAMPLE_SEP = ImmutableSampleData.of("sep-sample", SEP2015, null, 0D);
    private final static SampleData LIMS_SAMPLE_NOV = ImmutableSampleData.of("nov-sample", NOV2015, null, 0D);

    private final static SampleData LIMS_ARRIVED_NOV_SAMPLED_MAR = ImmutableSampleData.of("mar-sample-arrived-nov", NOV2015, MAR2015, 0D);

    private final static BiopsyData ECRF_BIOPSY_JAN = ImmutableBiopsyData.of(JAN2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData ECRF_BIOPSY_FEB = ImmutableBiopsyData.of(FEB2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData ECRF_BIOPSY_MAR = ImmutableBiopsyData.of(MAR2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData ECRF_BIOPSY_JUL = ImmutableBiopsyData.of(JUL2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData ECRF_BIOPSY_SEP = ImmutableBiopsyData.of(SEP2015, "", "", "", "", FormStatusState.UNKNOWN, false);

    private final static BiopsyData BIOPSY_NULL = ImmutableBiopsyData.of(null, "", "", "", "", FormStatusState.UNKNOWN, false);

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

    private static void runMatcherAndVerify(@NotNull List<SampleData> sequencedBiopsies, @NotNull List<BiopsyData> clinicalBiopsies,
            @Nullable String... expectedSampleIds) {
        final List<BiopsyData> matchedBiopsies = BiopsyMatcher.matchBiopsiesToTumorSamples("patient", sequencedBiopsies, clinicalBiopsies).values();
        matchedBiopsies.sort(CLINICAL_DATA_COMPARATOR);
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

    @NotNull
    private static final Comparator<BiopsyData> CLINICAL_DATA_COMPARATOR = (o1, o2) -> {
        final LocalDate o1Date = o1.date();
        final LocalDate o2Date = o2.date();
        if (o1Date == null) {
            return -1;
        } else if (o2Date == null) {
            return 1;
        } else {
            return o1Date.compareTo(o2Date);
        }
    };
}
