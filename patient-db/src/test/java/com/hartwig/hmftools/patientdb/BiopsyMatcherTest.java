package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
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

    private final static BiopsyClinicalData biopsyJan = new BiopsyClinicalData(jan2015, "", null);
    private final static BiopsyClinicalData biopsyFeb = new BiopsyClinicalData(feb2015, "", null);
    private final static BiopsyClinicalData biopsyMar = new BiopsyClinicalData(mar2015, "", null);
    private final static BiopsyClinicalData biopsyJul = new BiopsyClinicalData(jul2015, "", null);

    private final static BiopsyClinicalData biopsySep = new BiopsyClinicalData(sep2015, "", null);

    private final static BiopsyClinicalData biopsyNull = new BiopsyClinicalData(null, "", null);

    // MIVO:    ---biopsy(v)/sample---
    @Test
    public void testBiopsyAndSampleSameDate() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJul);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(v)---sample---
    @Test
    public void test1BiopsyBefore1SampleWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(x)---sample---
    @Test
    public void test1BiopsyBefore1SampleOutsideThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---sample---biopsy---
    @Test
    public void test1BiopsiesAfter1Sample() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsySep);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(v)---biopsy(v)---sample---
    @Test
    public void test2BiopsiesBefore1Sample2WithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb, biopsyMar);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(v)---null---sample---
    @Test
    public void test1Biopsy1NullBefore1SampleWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb, biopsyNull);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---null---biopsy(v)---sample---
    @Test
    public void test1Null1BiopsyBefore1SampleWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyNull, biopsyFeb);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(x)---biopsy(v)---sample---
    @Test
    public void test2BiopsiesBefore1Sample1WithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan, biopsyMar);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals("jul", matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(v)---biopsy(x)---sample---
    @Test
    public void test2BiopsiesBefore1Sample1WithinThresholdOutOfOrder() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar, biopsyJan);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(x)---biopsy(x)---sample---
    @Test
    public void test2BiopsiesBefore1Sample2OutsideThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsySep);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan, biopsyFeb);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(v)---sample---sample---
    @Test
    public void test1BiopsyBefore2SamplesWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul, limsBiopsyAug);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
    }

    // MIVO:    ---biopsy(v)---sample---biopsy---
    @Test
    public void test1SampleBetween2BiopsiesWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyJul);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar, biopsySep);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("jul", matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(v)---biopsy(v)---sample---sample---
    @Test
    public void test2SamplesAfter2Biopsies() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug, limsBiopsySep);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyFeb, biopsyMar);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(v)---sample---biopsy(v)---sample---
    @Test
    public void test2SamplesBetween2BiopsiesWithinThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug, limsBiopsyNov);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyMar, biopsySep);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        System.out.println(matchedBiopsies.size());
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals("aug", matchedBiopsies.get(0).sampleId());
        assertEquals("nov", matchedBiopsies.get(1).sampleId());
    }

    // MIVO:    ---biopsy(x)---sample---biopsy(v)---sample---
    @Test
    public void test2SamplesBetween2Biopsies1OutsideThreshold() {
        final List<BiopsyLimsData> sequencedBiopsies = Lists.newArrayList(limsBiopsyAug, limsBiopsyNov);
        final List<BiopsyClinicalData> clinicalBiopsies = Lists.newArrayList(biopsyJan, biopsySep);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies("patient", sequencedBiopsies,
                clinicalBiopsies);
        System.out.println(matchedBiopsies.size());
        assertTrue(clinicalBiopsies.size() == matchedBiopsies.size());
        assertEquals(null, matchedBiopsies.get(0).sampleId());
        assertEquals(null, matchedBiopsies.get(1).sampleId());
    }
}
