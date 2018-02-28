package com.hartwig.hmftools.patientdb.matchers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentDrugData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentMatcherTest {
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private final static LocalDate MAY2015 = LocalDate.parse("2015-05-01");
    private final static LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private final static LocalDate SEP2015 = LocalDate.parse("2015-09-01");

    private final static BiopsyTreatmentData TREATMENT_FEB_JUL2015 = ImmutableBiopsyTreatmentData.of("Yes",
            Lists.newArrayList(drugWithStartAndEndDate(FEB2015, JUL2015)),
            FormStatusState.UNKNOWN,
            false);
    private final static BiopsyTreatmentData TREATMENT_MAY_SEP2015 = ImmutableBiopsyTreatmentData.of("Yes",
            Lists.newArrayList(drugWithStartAndEndDate(MAY2015, SEP2015)),
            FormStatusState.UNKNOWN,
            false);
    private final static BiopsyTreatmentData NO_TREATMENT_GIVEN =
            ImmutableBiopsyTreatmentData.of("No", Lists.newArrayList(), FormStatusState.UNKNOWN, false);
    private final static BiopsyTreatmentData TREATMENT_MAR_NULL = ImmutableBiopsyTreatmentData.of("Yes",
            Lists.newArrayList(drugWithStartAndEndDate(MAR2015, null)),
            FormStatusState.UNKNOWN,
            false);

    private final static BiopsyData BIOPSY_JAN = ImmutableBiopsyData.of(JAN2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData BIOPSY_FEB = ImmutableBiopsyData.of(FEB2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData BIOPSY_MAR = ImmutableBiopsyData.of(MAR2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData BIOPSY_SEP = ImmutableBiopsyData.of(SEP2015, "", "", "", "", FormStatusState.UNKNOWN, false);
    private final static BiopsyData BIOPSY_NULL = ImmutableBiopsyData.of(null, "", "", "", "", FormStatusState.UNKNOWN, false);

    // MIVO:    ---start(feb)-biopsy(mar)----end(jul)---
    @Test
    public void testTreatmentStartBeforeBiopsyFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    // MIVO:    ---start/biopsy(feb)-----end(jul)---
    @Test
    public void testTreatmentStartSameDateBiopsySucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_FEB);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)-start(feb)-----end(jul)---
    @Test
    public void testTreatmentStartAfterBiopsySucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)----start(may)----end(sep)---
    @Test
    public void testTreatmentStart4MonthsAfterBiopsyFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAY_SEP2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    // MIVO:    ---biopsy(jan)-start(feb)-----end(jul)--biopsy(sep)-
    @Test
    public void testDoesntMatchTreatmentsWithTreatmentGivenNo() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015, NO_TREATMENT_GIVEN);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
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
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_FEB);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    @Test
    public void testDoesntMatchBiopsyWithNullDate() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_NULL);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    @NotNull
    private static BiopsyTreatmentDrugData drugWithStartAndEndDate(@Nullable LocalDate startDate, @Nullable LocalDate endDate) {
        return ImmutableBiopsyTreatmentDrugData.of("anything", startDate, endDate, Lists.newArrayList());
    }
}
