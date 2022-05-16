package com.hartwig.hmftools.patientdb.clinical.matchers;

import static com.hartwig.hmftools.patientdb.clinical.datamodel.ClinicalTestFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.clinical.datamodel.ClinicalTestFactory.biopsyTreatmentResponseBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.DrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableDrugData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentResponseMatcherTest {

    private static final LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private static final LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private static final LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private static final LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private static final LocalDate AUG2015 = LocalDate.parse("2015-08-01");
    private static final LocalDate SEP2015 = LocalDate.parse("2015-09-01");
    private static final LocalDate OCT2015 = LocalDate.parse("2015-10-01");
    private static final LocalDate NOV2015 = LocalDate.parse("2015-11-01");

    private static final BiopsyTreatmentData TREATMENT_JAN2015_JUL2015 = treatmentWithStartEnd(JAN2015, JUL2015);
    private static final BiopsyTreatmentData TREATMENT_FEB2015_JUL2015 = treatmentWithStartEnd(FEB2015, JUL2015);
    private static final BiopsyTreatmentData TREATMENT_FEB2015_NULL = treatmentWithStartEnd(FEB2015, null);
    private static final BiopsyTreatmentData TREATMENT_OCT2015_NULL = treatmentWithStartEnd(OCT2015, null);

    private static final BiopsyTreatmentResponseData RESPONSE_JAN2015 = responseOnDate(JAN2015);
    private static final BiopsyTreatmentResponseData RESPONSE_FEB2015 = responseOnDate(FEB2015);
    private static final BiopsyTreatmentResponseData RESPONSE_MAR2015 = responseOnDate(MAR2015);
    private static final BiopsyTreatmentResponseData RESPONSE_JUL2015 = responseOnDate(JUL2015);
    private static final BiopsyTreatmentResponseData RESPONSE_AUG2015 = responseOnDate(AUG2015);
    private static final BiopsyTreatmentResponseData RESPONSE_AUG2015_BASELINE = baselineResponseOnDate(AUG2015);
    private static final BiopsyTreatmentResponseData RESPONSE_SEP2015 = responseOnDate(SEP2015);
    private static final BiopsyTreatmentResponseData RESPONSE_NOV2015 = responseOnDate(NOV2015);
    private static final BiopsyTreatmentResponseData RESPONSE_NULL = responseOnDate(null);

    //    ---response(jan)-start(feb)-----end(jul)--
    @Test
    public void responseBeforeStartFails() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_JAN2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertNull(matchedResponses.get(0).treatmentId());
    }

    //    ---start/response(feb)-----end(jul)--
    @Test
    public void responseSameDateAsStartFails() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_FEB2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertNull(matchedResponses.get(0).treatmentId());
    }

    //    ---start(feb)-response(mar)----end(jul)--
    @Test
    public void responseAfterStartBeforeEndSucceeds() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    //    ---start(feb)-----end/response(jul)--
    @Test
    public void responseAfterStartSameDateAsEndSucceeds() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_JUL2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    //    ---start(feb)-----end(jul)-response(aug)--
    @Test
    public void responseAfterEndSucceeds() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_AUG2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    //    ---start(feb)-response(mar)------end(null)
    @Test
    public void responseAfterStartBeforeNullSucceeds() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    //    ---start(feb)-------end/response(null)
    @Test
    public void responseNullFails() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_NULL);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertNull(matchedResponses.get(0).treatmentId());
    }

    //    ---start(feb)-response(mar)----response(jul)--end(null)
    @Test
    public void multipleResponsesMatchTreatment() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015, RESPONSE_JUL2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());

        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
        assertMatch(treatments.get(0).id(), matchedResponses.get(1).treatmentId());
    }

    //    ---start1/start2(feb)-response(mar)----end1(jul)--end2(null)
    @Test
    public void overlappingTreatmentsCannotBeUsedForMatching() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015, TREATMENT_FEB2015_NULL);
        List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());
        assertNull(matchedResponses.get(0).treatmentId());
    }

    //    --start(jan)-response(feb)-response(mar)-end(jul)-responseNE(aug)-response(sep)-start(okt)-response(nov)
    @Test
    public void realTimelineWorksAsExpected() {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_JAN2015_JUL2015, TREATMENT_OCT2015_NULL);
        List<BiopsyTreatmentResponseData> responses =
                Lists.newArrayList(RESPONSE_FEB2015, RESPONSE_MAR2015, RESPONSE_AUG2015_BASELINE, RESPONSE_SEP2015, RESPONSE_NOV2015);
        List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertEquals(matchedResponses.size(), responses.size());

        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
        assertMatch(treatments.get(0).id(), matchedResponses.get(1).treatmentId());
        assertNull(matchedResponses.get(2).treatmentId());
        assertNull(matchedResponses.get(3).treatmentId());
        assertMatch(treatments.get(1).id(), matchedResponses.get(4).treatmentId());
    }

    private static void assertMatch(int expectedId, @Nullable Integer actualId) {
        assertNotNull(actualId);
        assertEquals(expectedId, actualId.intValue());
    }

    @NotNull
    private static BiopsyTreatmentData treatmentWithStartEnd(@Nullable LocalDate start, @Nullable LocalDate end) {
        DrugData drug = ImmutableDrugData.of("drug", start, end, null, Lists.newArrayList());
        return biopsyTreatmentBuilder().treatmentGiven("Yes").radiotherapyGiven("Yes").addDrugs(drug).build();
    }

    @NotNull
    private static BiopsyTreatmentResponseData responseOnDate(@Nullable LocalDate date) {
        return biopsyTreatmentResponseBuilder().responseDate(date).response("response").measurementDone("Yes").build();
    }

    @NotNull
    private static BiopsyTreatmentResponseData baselineResponseOnDate(@Nullable LocalDate date) {
        return biopsyTreatmentResponseBuilder().responseDate(date).response("NE").measurementDone("Yes").build();
    }
}
