package com.hartwig.hmftools.patientdb.matchers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;
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

    private final static BiopsyTreatmentData TREATMENT_FEB2015_JUL2015 = treatmentWithStartEnd(FEB2015, JUL2015);
    private final static BiopsyTreatmentData TREATMENT_FEB2015_NULL = treatmentWithStartEnd(FEB2015, null);

    private final static BiopsyTreatmentResponseData RESPONSE_JAN2015 = responseOnDate(JAN2015);
    private final static BiopsyTreatmentResponseData RESPONSE_FEB2015 = responseOnDate(FEB2015);
    private final static BiopsyTreatmentResponseData RESPONSE_MAR2015 = responseOnDate(MAR2015);
    private final static BiopsyTreatmentResponseData RESPONSE_JUL2015 = responseOnDate(JUL2015);
    private final static BiopsyTreatmentResponseData RESPONSE_AUG2015 = responseOnDate(AUG2015);
    private final static BiopsyTreatmentResponseData RESPONSE_NULL = responseOnDate(null);

    // MIVO:    ---response(jan)-start(feb)-----end(jul)--
    @Test
    public void testResponseBeforeStartFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_JAN2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start/response(feb)-----end(jul)--
    @Test
    public void testResponseSameDateAsStartFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_FEB2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-response(mar)----end(jul)--
    @Test
    public void testResponseAfterStartBeforeEndSucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        final Integer treatmentId = matchedResponses.get(0).treatmentId();
        assertNotNull(treatmentId);
        assertEquals(treatments.get(0).id(), treatmentId.intValue());
    }

    // MIVO:    ---start(feb)-----end/response(jul)--
    @Test
    public void testResponseAfterStartSameDateAsEndSucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_JUL2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        final Integer treatmentId = matchedResponses.get(0).treatmentId();
        assertNotNull(treatmentId);
        assertEquals(treatments.get(0).id(), treatmentId.intValue());
    }

    // MIVO:    ---start(feb)-----end(jul)-response(aug)--
    @Test
    public void testResponseAfterEndFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_AUG2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-response(mar)------end(null)
    @Test
    public void testResponseAfterStartBeforeNullSucceeds() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        final Integer treatmentId = matchedResponses.get(0).treatmentId();
        assertNotNull(treatmentId);
        assertEquals(treatments.get(0).id(), treatmentId.intValue());
    }

    // MIVO:    ---start(feb)-------end/response(null)
    @Test
    public void testResponseNullFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_NULL);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }

    // MIVO:    ---start(feb)-response(mar)----response(jul)--end(null)
    @Test
    public void testMultipleResponsesMatchTreatment() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015, RESPONSE_JUL2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        final Integer treatmentId1 = matchedResponses.get(0).treatmentId();
        final Integer treatmentId2 = matchedResponses.get(1).treatmentId();
        assertNotNull(treatmentId1);
        assertNotNull(treatmentId2);
        assertEquals(treatments.get(0).id(), treatmentId1.intValue());
        assertEquals(treatments.get(0).id(), treatmentId1.intValue());
    }

    // MIVO:    ---start1/start2(feb)-response(mar)----end1(jul)--end2(null)
    @Test
    public void testResponseMatchesMultipleFails() {
        final List<BiopsyTreatmentData> treatments = Lists.newArrayList(TREATMENT_FEB2015_JUL2015, TREATMENT_FEB2015_NULL);
        final List<BiopsyTreatmentResponseData> responses = Lists.newArrayList(RESPONSE_MAR2015);
        final List<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments("patient", treatments, responses).values();
        assertTrue(matchedResponses.size() == responses.size());
        assertEquals(null, matchedResponses.get(0).treatmentId());
    }


    @NotNull
    private static BiopsyTreatmentData treatmentWithStartEnd(@Nullable LocalDate start, @Nullable LocalDate end) {
        return ImmutableBiopsyTreatmentData.of("Yes", start, end, Lists.newArrayList(), FormStatusState.UNKNOWN, false);
    }

    @NotNull
    private static BiopsyTreatmentResponseData responseOnDate(@Nullable LocalDate date) {
        return ImmutableBiopsyTreatmentResponseData.of(date, date, "NE", "Yes", "No", FormStatusState.UNKNOWN, false);
    }
}
