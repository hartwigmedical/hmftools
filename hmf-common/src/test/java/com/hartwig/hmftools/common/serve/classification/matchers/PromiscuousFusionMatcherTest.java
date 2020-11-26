package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class PromiscuousFusionMatcherTest {

    @Test
    public void canAssessWhetherEventIsPromiscuousFusion() {
        EventMatcher matcher = new PromiscuousFusionMatcher();

        assertTrue(matcher.matches("BRAF", "BRAF fusion"));
        assertTrue(matcher.matches("ROS1", "ROS1 REARRANGEMENT"));

        assertFalse(matcher.matches("BRAF", "V600E"));
        assertFalse(matcher.matches("ALK", "EML4-ALK"));
    }
}