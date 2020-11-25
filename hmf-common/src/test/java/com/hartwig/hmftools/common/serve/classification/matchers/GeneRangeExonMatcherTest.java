package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class GeneRangeExonMatcherTest {

    @Test
    public void canAssessWhetherEventIsGeneRangeExon() {
        EventMatcher matcher = new GeneRangeExonMatcher();

        assertTrue(matcher.matches("EGFR", "RARE EX 18-21 MUT"));
        assertTrue(matcher.matches("EGFR", "EGFR exon 19 deletions"));
        assertTrue(matcher.matches("AXSL1", "EXON 12 MUTATION"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}