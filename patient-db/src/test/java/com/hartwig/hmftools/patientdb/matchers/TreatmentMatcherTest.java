package com.hartwig.hmftools.patientdb.matchers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;

import org.junit.Test;

public class TreatmentMatcherTest {
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private final static LocalDate MAY2015 = LocalDate.parse("2015-05-01");
    private final static LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private final static LocalDate SEP2015 = LocalDate.parse("2015-09-01");

    private final static BiopsyTreatmentData TREATMENT_FEB_JUL2015 = new BiopsyTreatmentData("Yes", FEB2015, JUL2015,
            Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_MAY_SEP2015 = new BiopsyTreatmentData("Yes", MAY2015, SEP2015,
            Lists.newArrayList());
    private final static BiopsyTreatmentData NO_TREATMENT_GIVEN = new BiopsyTreatmentData("No", null, null,
            Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_MAR_NULL = new BiopsyTreatmentData("Yes", MAR2015, null,
            Lists.newArrayList());

    private final static BiopsyClinicalData BIOPSY_JAN = new BiopsyClinicalData(JAN2015, "");
    private final static BiopsyClinicalData BIOPSY_FEB = new BiopsyClinicalData(FEB2015, "");
    private final static BiopsyClinicalData BIOPSY_MAR = new BiopsyClinicalData(MAR2015, "");
    private final static BiopsyClinicalData BIOPSY_SEP = new BiopsyClinicalData(SEP2015, "");
    private final static BiopsyClinicalData BIOPSY_NULL = new BiopsyClinicalData(null, "");

    // MIVO:    ---start(feb)-biopsy(mar)----end(jul)---
    @Test
    public void testTreatmentStartBeforeBiopsyFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies,
                treatments);
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    // MIVO:    ---start/biopsy(feb)-----end(jul)---
    @Test
    public void testTreatmentStartSameDateBiopsySucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList(BIOPSY_FEB);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies,
                treatments);
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)-start(feb)-----end(jul)---
    @Test
    public void testTreatmentStartAfterBiopsySucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies,
                treatments);
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)----start(may)----end(sep)---
    @Test
    public void testTreatmentStart4MonthsAfterBiopsyFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAY_SEP2015);
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies,
                treatments);
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    // MIVO:    ---biopsy(jan)-start(feb)-----end(jul)--biopsy(sep)-
    @Test
    public void testDoesntMatchTreatmentsWithTreatmentGivenNo() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015, NO_TREATMENT_GIVEN);
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies,
                treatments);
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId1.intValue());
        assertEquals(null, matchedTreatments.get(1).biopsyId());
    }

    // MIVO:    ---biopsy(jan)-biopsy(feb)-start(mar)-------end(null)
    @Test
    public void testDoesntMatchMultipleTreatments() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_FEB);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies,
                treatments);
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    @Test
    public void testDoesntMatchBiopsyWithNullDate() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        final List<BiopsyClinicalData> biopsies = Lists.newArrayList(BIOPSY_NULL);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies,
                treatments);
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

}
