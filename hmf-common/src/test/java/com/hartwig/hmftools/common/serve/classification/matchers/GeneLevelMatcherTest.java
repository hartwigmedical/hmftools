package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class GeneLevelMatcherTest {

    @Test
    public void canAssessWhetherEventIsGeneLevelEvent() {
        EventMatcher matcher = new GeneLevelMatcher();

        assertTrue(matcher.matches("AKT1", "AKT1 act mut"));
        assertTrue(matcher.matches("AKT1", "AKT1"));
        assertTrue(matcher.matches("ATK1", "positive"));
        assertTrue(matcher.matches("TP53", "biallelic inactivation"));

        assertFalse(matcher.matches("NPM", "Exon 12 mutation"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}