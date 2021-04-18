package com.hartwig.hmftools.patientreporter.germline;

import static org.junit.Assert.*;

import org.junit.Test;

public class GermlineConditionTest {

    @Test
    public void canExtractLimsAnalysisType() {
        assertEquals(GermlineCondition.NEVER, GermlineCondition.extractGermlineCondition("NEVER"));
        assertEquals(GermlineCondition.ALWAYS, GermlineCondition.extractGermlineCondition("ALWAYS"));
        assertEquals(GermlineCondition.ONLY_SPECIFIC_VARIANT, GermlineCondition.extractGermlineCondition("ONLY_SPECIFIC_VARIANT"));
        assertEquals(GermlineCondition.ONLY_GERMLINE_HOM, GermlineCondition.extractGermlineCondition("ONLY_GERMLINE_HOM"));

    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownCohortType() {
        //noinspection ResultOfMethodCallIgnored
        GermlineCondition.extractGermlineCondition("always");
    }

}