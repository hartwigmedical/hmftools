package com.hartwig.hmftools.patientreporter.germline;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GermlineConditionTest {

    @Test
    public void canExtractGermlineCondition() {
        assertEquals(GermlineCondition.NEVER, GermlineCondition.toGermlineCondition("NEVER"));
        assertEquals(GermlineCondition.ALWAYS, GermlineCondition.toGermlineCondition("ALWAYS"));
        assertEquals(GermlineCondition.ONLY_SPECIFIC_VARIANT, GermlineCondition.toGermlineCondition("ONLY_SPECIFIC_VARIANT"));
        assertEquals(GermlineCondition.ONLY_GERMLINE_HOM, GermlineCondition.toGermlineCondition("ONLY_GERMLINE_HOM"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownGermlineCondition() {
        GermlineCondition.toGermlineCondition("always");
    }
}