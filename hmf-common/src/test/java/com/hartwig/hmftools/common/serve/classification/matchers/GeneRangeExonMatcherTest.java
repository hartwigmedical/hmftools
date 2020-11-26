package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class GeneRangeExonMatcherTest {

    private static final String EXON_KEYWORD = "exon";
    private static final Set<String> EXON_RANGE_EVENTS = Sets.newHashSet("RARE EX 18-21 MUT");
    private static final Set<String> EXON_RANGE_KEYWORDS = Sets.newHashSet("deletion", "mutation");

    @Test
    public void canAssessWhetherEventIsGeneRangeExon() {
        EventMatcher matcher = new GeneRangeExonMatcher(EXON_KEYWORD, EXON_RANGE_EVENTS, EXON_RANGE_KEYWORDS);

        assertTrue(matcher.matches("EGFR", "RARE EX 18-21 MUT"));
        assertTrue(matcher.matches("EGFR", "EGFR exon 19 deletions"));
        assertTrue(matcher.matches("AXSL1", "EXON 12 MUTATION"));

        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}