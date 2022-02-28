package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;


public class HlaMatcherTest {

    private static final Set<String> HLA_EVENTS = Sets.newHashSet("hla");

    @Test
    public void canAssessWhetherEventIsHLA() {
        EventMatcher matcher = new HlaMatcher(HLA_EVENTS);

        assertTrue(matcher.matches("-", "hla"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }

}