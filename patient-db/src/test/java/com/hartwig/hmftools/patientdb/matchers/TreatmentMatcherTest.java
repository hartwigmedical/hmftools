package com.hartwig.hmftools.patientdb.matchers;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;

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

    private final static BiopsyData BIOPSY_JAN = biopsyBuilder().date(JAN2015).build();
    private final static BiopsyData BIOPSY_FEB = biopsyBuilder().date(FEB2015).build();
    private final static BiopsyData BIOPSY_MAR = biopsyBuilder().date(MAR2015).build();
    private final static BiopsyData BIOPSY_SEP = biopsyBuilder().date(SEP2015).build();
    private final static BiopsyData BIOPSY_NULL = biopsyBuilder().date(null).build();

    // LISC:    ---biopsy(mar)----no-treatment---
    @Test
    public void oneBiopsyNoTreatmentMatches() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(NO_TREATMENT_GIVEN);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---start(feb)-biopsy(mar)----end(jul)---
    @Test
    public void treatmentStartBeforeBiopsyFails() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        assertNull(matchedTreatments.get(0).biopsyId());
    }

    // MIVO:    ---start/biopsy(feb)-----end(jul)---
    @Test
    public void treatmentStartSameDateBiopsySucceeds() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_FEB);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)-start(feb)-----end(jul)---
    @Test
    public void treatmentStartAfterBiopsySucceeds() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        final Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    // MIVO:    ---biopsy(jan)----start(may)----end(sep)---
    @Test
    public void testTreatmentStart4MonthsAfterBiopsyFails() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAY_SEP2015);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        assertEquals(null, matchedTreatments.get(0).biopsyId());
    }

    // LISC:    ---biopsy(jan)-start(feb)--- end (jul) --- biopt(sep) --- no treatment
    @Test
    public void twoBiopsyMatchToTreatmentAndNoTreatment() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015, NO_TREATMENT_GIVEN);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        final Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        final Integer matchedBiopsyId2 = matchedTreatments.get(1).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertNotNull(matchedBiopsyId2);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId1.intValue());
        assertEquals(biopsies.get(1).id(), matchedBiopsyId2.intValue());
    }

    // LISC:    ---biopsy(jan)-no treatment --- start(feb)---end (jul) --- biopsy(sep)
    @Test
    public void twoBiopsyMatchToNoTreatmentAndTreatment() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(NO_TREATMENT_GIVEN, TREATMENT_FEB_JUL2015);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        final Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        final Integer matchedBiopsyId2 = matchedTreatments.get(1).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertNotNull(matchedBiopsyId2);
        assertEquals(biopsies.get(1).id(), matchedBiopsyId1.intValue());
        assertEquals(biopsies.get(0).id(), matchedBiopsyId2.intValue());
    }

    // MIVO:    ---biopsy(jan)-biopsy(feb)-start(mar)-------end(null)
    @Test
    public void matchesToMostRecentBiopsy() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_FEB);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        final Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertEquals(biopsies.get(1).id(), matchedBiopsyId1.intValue());
    }

    @Test
    public void doesntMatchBiopsyWithNullDate() {
        final List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_NULL);
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        final List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        assertNull(matchedTreatments.get(0).biopsyId());
    }

    @NotNull
    private static DrugData drugWithStartAndEndDate(@Nullable LocalDate startDate, @Nullable LocalDate endDate) {
        return ImmutableDrugData.of("anything", startDate, endDate, null, Lists.newArrayList());
    }

    @NotNull
    private static List<BiopsyTreatmentData> assertedMatch(@NotNull List<BiopsyData> biopsies,
            @NotNull List<BiopsyTreatmentData> treatments) {
        final List<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertTrue(treatments.size() == matchedTreatments.size());
        return matchedTreatments;
    }
}
