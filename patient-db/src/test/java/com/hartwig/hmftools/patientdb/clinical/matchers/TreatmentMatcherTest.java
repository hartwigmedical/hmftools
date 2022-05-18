package com.hartwig.hmftools.patientdb.clinical.matchers;

import static com.hartwig.hmftools.patientdb.clinical.datamodel.ClinicalTestFactory.biopsyBuilder;
import static com.hartwig.hmftools.patientdb.clinical.datamodel.ClinicalTestFactory.biopsyTreatmentBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.DrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableDrugData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentMatcherTest {

    private static final LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private static final LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private static final LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private static final LocalDate MAY2015 = LocalDate.parse("2015-05-01");
    private static final LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private static final LocalDate SEP2015 = LocalDate.parse("2015-09-01");

    private static final BiopsyTreatmentData TREATMENT_FEB_JUL2015 =
            biopsyTreatmentBuilder().treatmentGiven("Yes").addDrugs(drugWithStartAndEndDate(FEB2015, JUL2015)).build();
    private static final BiopsyTreatmentData TREATMENT_MAY_SEP2015 =
            biopsyTreatmentBuilder().treatmentGiven("Yes").addDrugs(drugWithStartAndEndDate(MAY2015, SEP2015)).build();
    private static final BiopsyTreatmentData NO_TREATMENT_GIVEN = biopsyTreatmentBuilder().treatmentGiven("No").build();
    private static final BiopsyTreatmentData TREATMENT_MAR_NULL =
            biopsyTreatmentBuilder().treatmentGiven("Yes").addDrugs(drugWithStartAndEndDate(MAR2015, null)).build();

    private static final BiopsyData BIOPSY_JAN = biopsyBuilder().date(JAN2015).build();
    private static final BiopsyData BIOPSY_FEB = biopsyBuilder().date(FEB2015).build();
    private static final BiopsyData BIOPSY_MAR = biopsyBuilder().date(MAR2015).build();
    private static final BiopsyData BIOPSY_SEP = biopsyBuilder().date(SEP2015).build();
    private static final BiopsyData BIOPSY_NULL = biopsyBuilder().date(null).build();

    //    ---biopsy(mar)----no-treatment---
    @Test
    public void oneBiopsyNoTreatmentMatches() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(NO_TREATMENT_GIVEN);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    //    ---start(feb)-biopsy(mar)----end(jul)---
    @Test
    public void treatmentStartBeforeBiopsyFails() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_MAR);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        assertNull(matchedTreatments.get(0).biopsyId());
    }

    //    ---start/biopsy(feb)-----end(jul)---
    @Test
    public void treatmentStartSameDateBiopsySucceeds() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_FEB);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    //    ---biopsy(jan)-start(feb)-----end(jul)---
    @Test
    public void treatmentStartAfterBiopsySucceeds() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        Integer matchedBiopsyId = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId.intValue());
    }

    //    ---biopsy(jan)----start(may)----end(sep)---
    @Test
    public void testTreatmentStart4MonthsAfterBiopsyFails() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAY_SEP2015);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        assertNull(matchedTreatments.get(0).biopsyId());
    }

    //    ---biopsy(jan)-start(feb)--- end (jul) --- biopsy(sep) --- no treatment
    @Test
    public void twoBiopsyMatchToTreatmentAndNoTreatment() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB_JUL2015, NO_TREATMENT_GIVEN);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        Integer matchedBiopsyId2 = matchedTreatments.get(1).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertNotNull(matchedBiopsyId2);
        assertEquals(biopsies.get(0).id(), matchedBiopsyId1.intValue());
        assertEquals(biopsies.get(1).id(), matchedBiopsyId2.intValue());
    }

    //    ---biopsy(jan)-no treatment --- start(feb)---end (jul) --- biopsy(sep)
    @Test
    public void twoBiopsyMatchToNoTreatmentAndTreatment() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_SEP);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(NO_TREATMENT_GIVEN, TREATMENT_FEB_JUL2015);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        Integer matchedBiopsyId2 = matchedTreatments.get(1).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertNotNull(matchedBiopsyId2);
        assertEquals(biopsies.get(1).id(), matchedBiopsyId1.intValue());
        assertEquals(biopsies.get(0).id(), matchedBiopsyId2.intValue());
    }

    //    ---biopsy(jan)-biopsy(feb)-start(mar)-------end(null)
    @Test
    public void matchesToMostRecentBiopsy() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_JAN, BIOPSY_FEB);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        Integer matchedBiopsyId1 = matchedTreatments.get(0).biopsyId();
        assertNotNull(matchedBiopsyId1);
        assertEquals(biopsies.get(1).id(), matchedBiopsyId1.intValue());
    }

    //    --- biopsy (null) - start(mar) ----
    @Test
    public void doesntMatchBiopsyWithNullDate() {
        List<BiopsyData> biopsies = Lists.newArrayList(BIOPSY_NULL);
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_MAR_NULL);
        List<BiopsyTreatmentData> matchedTreatments = assertedMatch(biopsies, treatments);

        assertNull(matchedTreatments.get(0).biopsyId());
    }

    @NotNull
    private static DrugData drugWithStartAndEndDate(@Nullable LocalDate startDate, @Nullable LocalDate endDate) {
        return ImmutableDrugData.of("anything", startDate, endDate, null, Lists.newArrayList());
    }

    @NotNull
    private static List<BiopsyTreatmentData> assertedMatch(@NotNull List<BiopsyData> biopsies,
            @NotNull List<BiopsyTreatmentData> treatments) {
        List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies("patient", biopsies, treatments).values();
        assertEquals(treatments.size(), matchedTreatments.size());
        return matchedTreatments;
    }
}
