package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SignatureMatcherTest {

    @Test
    public void canAssessWhetherEventIsSignature() {
        EventMatcher matcher = new SignatureMatcher();

        assertTrue(matcher.matches("-", "Microsatellite Instability-High"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}