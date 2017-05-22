package com.hartwig.hmftools.patientdb.matchers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BiopsyMatcherTest {
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private final static LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private final static LocalDate AUG2015 = LocalDate.parse("2015-08-01");
    private final static LocalDate SEP2015 = LocalDate.parse("2015-09-01");
    private final static LocalDate NOV2015 = LocalDate.parse("2015-11-01");

    private final static BiopsyLimsData LIMS_BIOPSY_JUL = new BiopsyLimsData("jul", JUL2015, null);
    private final static BiopsyLimsData LIMS_BIOPSY_AUG = new BiopsyLimsData("aug", AUG2015, null);
    private final static BiopsyLimsData LIMS_BIOPSY_SEP = new BiopsyLimsData("sep", SEP2015, null);
    private final static BiopsyLimsData LIMS_BIOPSY_NOV = new BiopsyLimsData("nov", NOV2015, null);

    private final static BiopsyClinicalData BIOPSY_JAN = new BiopsyClinicalData(JAN2015, "");
    private final static BiopsyClinicalData BIOPSY_FEB = new BiopsyClinicalData(FEB2015, "");
    private final static BiopsyClinicalData BIOPSY_MAR = new BiopsyClinicalData(MAR2015, "");
    private final static BiopsyClinicalData BIOPSY_JUL = new BiopsyClinicalData(JUL2015, "");

    private final static BiopsyClinicalData BIOPSY_SEP = new BiopsyClinicalData(SEP2015, "");

    private final static BiopsyClinicalData BIOPSY_NULL = new BiopsyClinicalData(null, "");

    // MIVO:    ---biopsy(jul)/sample(jul)---
    @Test
    public void testBiopsyAndSampleSameDate() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_JUL);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void test1BiopsyBefore1SampleWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_MAR);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-x-x-x-x-sample(jul)---
    @Test
    public void test1BiopsyBefore1SampleOutsideThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---sample(aug)-biopsy(sep)---
    @Test
    public void test1BiopsiesAfter1Sample() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_AUG);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_SEP);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(feb)-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void test2BiopsiesBefore1Sample2WithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_FEB, BIOPSY_MAR);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(feb)-x-x-null(?)-x-sample(jul)---
    @Test
    public void test1Biopsy1NullBefore1SampleWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_FEB, BIOPSY_NULL);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---null(?)-biopsy(feb)-x-x-x-x-sample(jul)---
    @Test
    public void test1Null1BiopsyBefore1SampleWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_NULL, BIOPSY_FEB);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void test2BiopsiesBefore1Sample1WithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_MAR);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals("jul", matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void test2BiopsiesBefore1Sample1WithinThresholdOutOfOrder() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_MAR, BIOPSY_JAN);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals("jul", matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(jan)-biopsy(feb)-x-x-x-x-x-x-sample(sep)---
    @Test
    public void test2BiopsiesBefore1Sample2OutsideThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_SEP);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_FEB);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(feb)-x-x-x-x-sample(jul)-sample(aug)---
    @Test
    public void test1BiopsyBefore2SamplesWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL, LIMS_BIOPSY_AUG);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_FEB);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(mar)-x-x-x-sample(jul)-x-biopsy(sep)---
    @Test
    public void test1SampleBetween2BiopsiesWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_JUL);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_MAR, BIOPSY_SEP);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(feb)-biopsy(mar)-x-x-x-x-sample(aug)-sample(sep)---
    @Test
    public void test2SamplesAfter2Biopsies() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_AUG, LIMS_BIOPSY_SEP);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_FEB, BIOPSY_MAR);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(mar)-x-x-x-x-sample(aug)-biopsy(sep)-x-sample(nov)---
    @Test
    public void test2SamplesBetween2BiopsiesWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_AUG, LIMS_BIOPSY_NOV);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_MAR, BIOPSY_SEP);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("aug", matchedBiopsies.get(0).sampleId());
        assertEquals("nov", matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-x-x-x-x-x-sample(aug)-biopsy(sep)-x-sample(nov)---
    @Test
    public void test2SamplesBetween2Biopsies1OutsideThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(LIMS_BIOPSY_AUG, LIMS_BIOPSY_NOV);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    @NotNull
    private static Comparator<BiopsyClinicalData> clinicalDataComparator = (o1, o2) -> {
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
