package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FusionPairAndExonRangeMatcherTest {

    @Test
    public void canAssessWhetherEventIsFusionAndExonRange() {
        EventMatcher matcher = new FusionPairAndExonRangeMatcher();

        assertTrue(matcher.matches("MET", "EXON 14 SKIPPING MUTATION"));
        assertFalse(matcher.matches("NRG1", "EXON 14 SKIPPING MUTATION"));
    }
}