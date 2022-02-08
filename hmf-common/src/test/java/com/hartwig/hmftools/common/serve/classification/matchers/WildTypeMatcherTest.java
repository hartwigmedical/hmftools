package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class WildTypeMatcherTest {

    private static final Set<String> WILD_TYPE_GENE_KEY_PHRASES = Sets.newHashSet("wild-type");

    @Test
    public void canAssessWhetherEventIsWildType() {
        WildTypeMatcher matcher = new WildTypeMatcher(WILD_TYPE_GENE_KEY_PHRASES);

        assertTrue(matcher.matches("EGFR", "EGFR wild-type"));

        assertFalse(matcher.matches("EGFR", "EGFR wild"));

    }

}