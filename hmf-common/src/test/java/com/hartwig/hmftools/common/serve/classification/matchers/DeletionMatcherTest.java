package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class DeletionMatcherTest {

    @Test
    public void canAssessWhetherEventIsDeletion() {
        EventMatcher matcher = new DeletionMatcher();

        assertTrue(matcher.matches("CDKN2A", "CDKN2A del"));
        assertTrue(matcher.matches("CDKN2A", "CDKN2A dec exp"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}