package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;
import com.hartwig.hmftools.patientdb.matchers.BiopsyMatcher;

import org.junit.Test;

public class BiopsyMatcherTest {
    private final static LocalDate jan2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate feb2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate mar2015 = LocalDate.parse("2015-03-01");
    private final static LocalDate jul2015 = LocalDate.parse("2015-07-01");
    private final static LocalDate aug2015 = LocalDate.parse("2015-08-01");
    private final static LocalDate sep2015 = LocalDate.parse("2015-09-01");
    private final static LocalDate nov2015 = LocalDate.parse("2015-11-01");

    private final static BiopsyLimsData limsBiopsyJul = new BiopsyLimsData("jul", jul2015, null);
    private final static BiopsyLimsData limsBiopsyAug = new BiopsyLimsData("aug", aug2015, null);
    private final static BiopsyLimsData limsBiopsySep = new BiopsyLimsData("sep", sep2015, null);
    private final static BiopsyLimsData limsBiopsyNov = new BiopsyLimsData("nov", nov2015, null);

    private final static BiopsyClinicalData biopsyJan = new BiopsyClinicalData(jan2015, "");
    private final static BiopsyClinicalData biopsyFeb = new BiopsyClinicalData(feb2015, "");
    private final static BiopsyClinicalData biopsyMar = new BiopsyClinicalData(mar2015, "");
    private final static BiopsyClinicalData biopsyJul = new BiopsyClinicalData(jul2015, "");

    private final static BiopsyClinicalData biopsySep = new BiopsyClinicalData(sep2015, "");

    private final static BiopsyClinicalData biopsyNull = new BiopsyClinicalData(null, "");

    // MIVO:    ---biopsy(jul)/sample(jul)---
    @Test
    public void testBiopsyAndSampleSameDate() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJul);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void test1BiopsyBefore1SampleWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(jan)-x-x-x-x-x-sample(jul)---
    @Test
    public void test1BiopsyBefore1SampleOutsideThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---sample(aug)-biopsy(sep)---
    @Test
    public void test1BiopsiesAfter1Sample() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsySep);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(feb)-biopsy(mar)-x-x-x-sample(jul)---
    @Test
    public void test2BiopsiesBefore1Sample2WithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb, biopsyMar);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb, biopsyNull);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyNull, biopsyFeb);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan, biopsyMar);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar, biopsyJan);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsySep);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan, biopsyFeb);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul, limsBiopsyAug);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(mar)-x-x-x-sample(jul)-x-biopsy(sep)---
    @Test
    public void test1SampleBetween2BiopsiesWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar, biopsySep);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug, limsBiopsySep);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb, biopsyMar);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug, limsBiopsyNov);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar, biopsySep);
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
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug, limsBiopsyNov);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan, biopsySep);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        matchedBiopsies.sort(clinicalDataComparator);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

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
