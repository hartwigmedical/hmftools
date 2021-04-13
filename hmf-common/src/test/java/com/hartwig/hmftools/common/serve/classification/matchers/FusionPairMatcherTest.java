package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class FusionPairMatcherTest {

    private static final Set<String> EXONIC_DEL_DUP_FUSION_KEY_PHRASES = Sets.newHashSet("exonic del dup");
    private static final Set<String> EXONIC_DEL_DUP_FUSION_EVENTS = Sets.newHashSet("VIII");
    private static final Set<String> FUSION_EVENTS_TO_SKIP = Sets.newHashSet("AR-V7");

    @Test
    public void canAssessWhetherEventIsFusionPair() {
        EventMatcher matcher =
                new FusionPairMatcher(EXONIC_DEL_DUP_FUSION_KEY_PHRASES, EXONIC_DEL_DUP_FUSION_EVENTS, FUSION_EVENTS_TO_SKIP);

        assertTrue(matcher.matches("ABL1", "ABL1-BCR fusion"));
        assertTrue(matcher.matches("ALK", "EML4-ALK"));
        assertTrue(matcher.matches("EGFR", "VIII"));
        assertTrue(matcher.matches("EGFR", "exonic del dup of any exon!"));
        assertTrue(matcher.matches("NKX2-1", "TRB-NKX2-1 Fusion"));

        assertFalse(matcher.matches("BRAF", "BRAF-"));
        assertFalse(matcher.matches("BRAF", "BRAF amp"));
        assertFalse(matcher.matches("AR", "AR-V7"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}