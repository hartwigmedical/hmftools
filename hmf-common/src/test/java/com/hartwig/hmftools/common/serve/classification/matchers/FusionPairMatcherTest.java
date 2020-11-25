package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class FusionPairMatcherTest {

    @Test
    public void canAssessWhetherEventIsFusionPair() {
        EventMatcher matcher = new FusionPairMatcher();

        assertTrue(matcher.matches("ABL1", "ABL1-BCR fusion"));
        assertTrue(matcher.matches("ALK", "EML4-ALK"));
        assertTrue(matcher.matches("EGFR", "VIII"));
        assertTrue(matcher.matches("NKX2-1", "TRB-NKX2-1 Fusion"));

        assertFalse(matcher.matches("BRAF", "BRAF amp"));
        assertFalse(matcher.matches("AR", "AR-V7"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}