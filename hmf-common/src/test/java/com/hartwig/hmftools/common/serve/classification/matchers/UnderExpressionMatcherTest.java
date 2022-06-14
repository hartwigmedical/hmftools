package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

import junit.framework.TestCase;

public class UnderExpressionMatcherTest {

    private static final Set<String> UNDER_EXPRESSION_KEYWORDS = Sets.newHashSet();
    private static final Set<String> UNDER_EXPRESSION_KEY_PHRASES = Sets.newHashSet("dec exp");

    @Test
    public void canAssessWhetherEventIsDeletion() {
        EventMatcher matcher = new UnderExpressionMatcher(
                UNDER_EXPRESSION_KEYWORDS,
                UNDER_EXPRESSION_KEY_PHRASES);

        assertTrue(matcher.matches("CDKN2A", "CDKN2A dec exp"));

        assertFalse(matcher.matches("CDKN2A", "CDKN2A del"));
        assertFalse(matcher.matches("EGFR", "EGFR Ex19 del"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }

}