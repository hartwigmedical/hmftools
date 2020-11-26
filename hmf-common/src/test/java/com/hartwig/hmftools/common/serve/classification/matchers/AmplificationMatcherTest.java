package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class AmplificationMatcherTest {

    @Test
    public void canAssessWhetherEventIsAmplification() {
        EventMatcher matcher = new AmplificationMatcher();

        assertTrue(matcher.matches("ALK", "ALK over exp"));
        assertTrue(matcher.matches("ALK", "ALK  amp"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}