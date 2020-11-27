package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class SignatureMatcherTest {

    private static final Set<String> SIGNATURE_EVENTS = Sets.newHashSet("Microsatellite Instability-High");

    @Test
    public void canAssessWhetherEventIsSignature() {
        EventMatcher matcher = new SignatureMatcher(SIGNATURE_EVENTS);

        assertTrue(matcher.matches("-", "Microsatellite Instability-High"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}