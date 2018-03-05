package com.hartwig.hmftools.patientdb.matchers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentResponseMatcherTest {
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");
    private final static LocalDate JUL2015 = LocalDate.parse("2015-07-01");
    private final static LocalDate AUG2015 = LocalDate.parse("2015-08-01");
    private final static LocalDate SEP2015 = LocalDate.parse("2015-09-01");
    private final static LocalDate OCT2015 = LocalDate.parse("2015-10-01");
    private final static LocalDate NOV2015 = LocalDate.parse("2015-11-01");

    private final static BiopsyTreatmentData TREATMENT_JAN2015_JUL2015 = treatmentWithStartEnd(JAN2015, JUL2015);
    private final static BiopsyTreatmentData TREATMENT_FEB2015_JUL2015 = treatmentWithStartEnd(FEB2015, JUL2015);
    private final static BiopsyTreatmentData TREATMENT_FEB2015_NULL = treatmentWithStartEnd(FEB2015, null);
    private final static BiopsyTreatmentData TREATMENT_OCT2015_NULL = treatmentWithStartEnd(OCT2015, null);

    private final static BiopsyTreatmentResponseData RESPONSE_JAN2015 = responseOnDate(JAN2015);
    private final static BiopsyTreatmentResponseData RESPONSE_FEB2015 = responseOnDate(FEB2015);
    private final static BiopsyTreatmentResponseData RESPONSE_MAR2015 = responseOnDate(MAR2015);
    private final static BiopsyTreatmentResponseData RESPONSE_JUL2015 = responseOnDate(JUL2015);
    private final static BiopsyTreatmentResponseData RESPONSE_AUG2015 = responseOnDate(AUG2015);
    private final static BiopsyTreatmentResponseData RESPONSE_AUG2015_BASELINE = baselineResponseOnDate(AUG2015);
    private final static BiopsyTreatmentResponseData RESPONSE_SEP2015 = responseOnDate(SEP2015);
    private final static BiopsyTreatmentResponseData RESPONSE_NOV2015 = responseOnDate(NOV2015);
    private final static BiopsyTreatmentResponseData RESPONSE_NULL = responseOnDate(null);

    // MIVO:    ---response(jan)-start(feb)-----end(jul)--
    @Test
    public void responseBeforeStartFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_JAN2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start/response(feb)-----end(jul)--
    @Test
    public void responseSameDateAsStartFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_FEB2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-response(mar)----end(jul)--
    @Test
    public void responseAfterStartBeforeEndSucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-----end/response(jul)--
    @Test
    public void responseAfterStartSameDateAsEndSucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_JUL2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-----end(jul)-response(aug)--
    @Test
    public void responseAfterEndSucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_AUG2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-response(mar)------end(null)
    @Test
    public void responseAfterStartBeforeNullSucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-------end/response(null)
    @Test
    public void responseNullFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_NULL);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-response(mar)----response(jul)--end(null)
    @Test
    public void multipleResponsesMatchTreatment() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015, RESPONSE_JUL2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());

        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
        assertMatch(treatments.get(0).id(), matchedResponses.get(1).treatmentId());
    }

    // MIVO:    ---start1/start2(feb)-response(mar)----end1(jul)--end2(null)
    @Test
    public void overlappingTreatmentsCannotBeUsedForMatching() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015, TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // KODU:    --start(jan)-response(feb)-response(mar)-end(jul)-responseNE(aug)-response(sep)-start(okt)-response(nov)
    @Test
    public void realTimelineWorksAsExpected() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_JAN2015_JUL2015, TREATMENT_OCT2015_NULL);
        final List<BiopsyTreatmentResponseData> responses =
                Lists.newArrayList(RESPONSE_FEB2015, RESPONSE_MAR2015, RESPONSE_AUG2015_BASELINE, RESPONSE_SEP2015, RESPONSE_NOV2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());

        assertMatch(treatments.get(0).id(), matchedResponses.get(0).treatmentId());
        assertMatch(treatments.get(0).id(), matchedResponses.get(1).treatmentId());
        assertEquals(null, matchedResponses.get(2).treatmentId());
        assertEquals(null, matchedResponses.get(3).treatmentId());
        assertMatch(treatments.get(1).id(), matchedResponses.get(4).treatmentId());
    }

    private static void assertMatch(int expectedId, @Nullable Integer actualId) {
        assertNotNull(actualId);
        assertEquals(expectedId, actualId.intValue());
    }

    @NotNull
    private static BiopsyTreatmentData treatmentWithStartEnd(@Nullable LocalDate start, @Nullable LocalDate end) {
        BiopsyTreatmentDrugData drug = ImmutableBiopsyTreatmentDrugData.of("drug", start, end, Lists.newArrayList());
        return ImmutableBiopsyTreatmentData.of("Yes", Lists.newArrayList(drug), FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private static BiopsyTreatmentResponseData responseOnDate(@Nullable LocalDate date) {
        return ImmutableBiopsyTreatmentResponseData.of(date, date, "response", "Yes", "No", FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private static BiopsyTreatmentResponseData baselineResponseOnDate(@Nullable LocalDate date) {
        return ImmutableBiopsyTreatmentResponseData.of(date, date, "NE", "Yes", "No", FormStatusState.UNKNOWN, false);
    }
}
