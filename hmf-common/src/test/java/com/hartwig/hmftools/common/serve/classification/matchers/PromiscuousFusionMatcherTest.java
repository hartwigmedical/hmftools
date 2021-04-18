package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class PromiscuousFusionMatcherTest {

    private static final Set<String> PROMISCUOUS_FUSION_KEY_PHRASES = Sets.newHashSet("fusion", "REARRANGEMENT");

    @Test
    public void canAssessWhetherEventIsPromiscuousFusion() {
        EventMatcher matcher = new PromiscuousFusionMatcher(PROMISCUOUS_FUSION_KEY_PHRASES,
                new FusionPairMatcher(Sets.newHashSet(), Sets.newHashSet(), Sets.newHashSet()));

        assertTrue(matcher.matches("BRAF", "BRAF fusion"));
        assertTrue(matcher.matches("ROS1", "ROS1 REARRANGEMENT"));

        assertFalse(matcher.matches("BRAF", "V600E"));
        assertFalse(matcher.matches("ALK", "EML4-ALK fusion"));
    }
}