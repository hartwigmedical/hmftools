package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class DeletionMatcherTest {

    private static final Set<String> BLACKLIST_KEY_PHRASES = Sets.newHashSet("Ex19");
    private static final Set<String> DELETION_KEYWORDS = Sets.newHashSet("del");
    private static final Set<String> DELETION_KEY_PHRASES = Sets.newHashSet();

    @Test
    public void canAssessWhetherEventIsDeletion() {
        EventMatcher matcher = new DeletionMatcher(BLACKLIST_KEY_PHRASES,
                DELETION_KEYWORDS,
                DELETION_KEY_PHRASES);

        assertTrue(matcher.matches("CDKN2A", "CDKN2A del"));

        assertFalse(matcher.matches("CDKN2A", "CDKN2A dec exp"));
        assertFalse(matcher.matches("EGFR", "EGFR Ex19 del"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}