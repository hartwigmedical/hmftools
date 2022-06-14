package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class AmplificationMatcherTest {

    private static final Set<String> AMPLIFICATION_KEYWORDS = Sets.newHashSet("amp");
    private static final Set<String> AMPLIFICATION_KEY_PHRASES = Sets.newHashSet();
    private static final Set<String> OVER_EXPRESSION_KEYWORDS = Sets.newHashSet();
    private static final Set<String> OVER_EXPRESSION_KEY_PHRASES = Sets.newHashSet("over exp");

    @Test
    public void canAssessWhetherEventIsAmplification() {
        EventMatcher matcher = new AmplificationMatcher(AMPLIFICATION_KEYWORDS,
                AMPLIFICATION_KEY_PHRASES,
                OVER_EXPRESSION_KEYWORDS,
                OVER_EXPRESSION_KEY_PHRASES);

        assertTrue(matcher.matches("ALK", "ALK over exp"));
        assertTrue(matcher.matches("ALK", "ALK  amp"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}