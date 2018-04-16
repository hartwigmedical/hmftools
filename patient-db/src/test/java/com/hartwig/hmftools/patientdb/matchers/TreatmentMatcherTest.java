package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.sampleBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.data.SampleData;

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

    private final static BiopsyTreatmentData TREATMENT_FEB_JUL2015 =
            biopsyTreatmentBuilder().treatmentGiven("Yes").addDrugs(drugWithStartAndEndDate(FEB2015, JUL2015)).build();
    private final static BiopsyTreatmentData TREATMENT_MAY_SEP2015 =
            biopsyTreatmentBuilder().treatmentGiven("Yes").addDrugs(drugWithStartAndEndDate(MAY2015, SEP2015)).build();
    private final static BiopsyTreatmentData NO_TREATMENT_GIVEN = biopsyTreatmentBuilder().treatmentGiven("No").build();
    private final static BiopsyTreatmentData TREATMENT_MAR_NULL =
            biopsyTreatmentBuilder().treatmentGiven("Yes").addDrugs(drugWithStartAndEndDate(MAR2015, null)).build();

    private final static SampleData LIMS_SAMPLE_JUL = sampleBuilder(JUL2015).build();

    private final static BiopsyData BIOPSY_JAN = biopsyBuilder().date(JAN2015).build();
    private final static BiopsyData BIOPSY_FEB = biopsyBuilder().date(FEB2015).build();
    private final static BiopsyData BIOPSY_MAR = biopsyBuilder().date(MAR2015).build();
    private final static BiopsyData BIOPSY_SEP = biopsyBuilder().date(SEP2015).build();
    private final static BiopsyData BIOPSY_NULL = biopsyBuilder().date(null).build();

    // LISC:    ---biopsy(mar)----no-treatment---
    @Test
    public void testOneBiopsyNoTreatment() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(NO_TREATMENT_GIVEN);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---start(feb)-biopsy(mar)----end(jul)---
    @Test
    public void testTreatmentStartBeforeBiopsyFails() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    // MIVO:    ---start/biopsy(feb)-----end(jul)---
    @Test
    public void testTreatmentStartSameDateBiopsySucceeds() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_FEB);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)-start(feb)-----end(jul)---
    @Test
    public void testTreatmentStartAfterBiopsySucceeds() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)----start(may)----end(sep)---
    @Test
    public void testTreatmentStart4MonthsAfterBiopsyFails() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAY_SEP2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    // LISC:    ---biopsy(jan)-start(feb)---end (jul) - biopt(sep) --- no treatment
    @Test
    public void testTwoBiopsyMatchToTreatmentAndNoTreatment() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015, NO_TREATMENT_GIVEN);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        final Integer matchedBiopsyId2 = matchedTreatments.get(1).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertNotNull(matchedBiopsyId2);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId1.intValue());
        assertEquals(biopsies.get(1).id(), matchedBiopsyId2.intValue());
    }

    // LISC:    ---biopsy(jan)-no treatment --- start(feb)---end (jul) --- biopsy(sep)
    @Test
    public void testTwoBiopsyMatchToNoTreatmentAndTreatment() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(NO_TREATMENT_GIVEN, TREATMENT_FEB_JUL2015);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        final Integer matchedBiopsyId2 = matchedTreatments.get(1).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertNotNull(matchedBiopsyId2);
        assertEquals(biopsies.get(1).id(), matchedBiopsyId1.intValue());
        assertEquals(biopsies.get(0).id(), matchedBiopsyId2.intValue());
    }

    // MIVO:    ---biopsy(jan)-start(feb)-----end(jul)--biopsy(sep)-
    @Test
    public void testDoesntMatchTreatmentsWithTreatmentGivenNo() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015, NO_TREATMENT_GIVEN);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        final Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId1.intValue());
        assertEquals(matchedBiopsyId1, matchedTreatments.get(1).biopsyId());
    }

    // MIVO:    ---biopsy(jan)-biopsy(feb)-start(mar)-------end(null)
    @Test
    public void testDoesntMatchMultipleTreatments() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_FEB);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    @Test
    public void testDoesntMatchBiopsyWithNullDate() {
        final List<SampleData> sequencedBiopsies = Lists.newArrayList(LIMS_SAMPLE_JUL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_NULL);
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments, sequencedBiopsies).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    @NotNull
    private static DrugData drugWithStartAndEndDate(@Nullable LocalDate startDate, @Nullable LocalDate endDate) {
        return ImmutableDrugData.of("anything", startDate, endDate, null, Lists.newArrayList());
    }
}
