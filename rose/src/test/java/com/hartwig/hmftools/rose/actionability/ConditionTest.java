package com.hartwig.hmftools.rose.actionability;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class ConditionTest {

    @Test
    public void canExtractCondition() {
        assertEquals(Condition.ONLY_HIGH, Condition.toCondition("ONLY_HIGH"));
        assertEquals(Condition.ALWAYS, Condition.toCondition("ALWAYS"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownCondition() {
        Condition.toCondition("always");
    }
}