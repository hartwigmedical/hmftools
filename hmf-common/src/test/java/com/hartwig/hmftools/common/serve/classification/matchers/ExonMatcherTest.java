package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class ExonMatcherTest {

    private static final Set<String> EXON_IDENTIFIERS = Sets.newHashSet("exon", "EXON");
    private static final Set<String> EXON_KEYWORDS = Sets.newHashSet("deletions", "MUTATION");
    private static final Set<String> EXON_BLACKLIST_KEY_PHRASES= Sets.newHashSet("no exon");
    private static final Set<String> SPECIFIC_EXON_EVENTS = Sets.newHashSet("RARE EX 18-21 MUT");

    @Test
    public void canAssessWhetherEventIsExon() {
        EventMatcher matcher = new ExonMatcher(EXON_IDENTIFIERS, EXON_KEYWORDS, EXON_BLACKLIST_KEY_PHRASES, SPECIFIC_EXON_EVENTS);

        assertTrue(matcher.matches("EGFR", "EGFR exon 19 deletions"));
        assertTrue(matcher.matches("AXSL1", "EXON 12 MUTATION"));
        assertTrue(matcher.matches("EGFR", "RARE EX 18-21 MUT"));

        assertFalse(matcher.matches("EGFR", "EGFR no exon deletions"));
        assertFalse(matcher.matches("EGFR", "EGFR deletions"));
        assertFalse(matcher.matches("EGFR", "EGFR exon 19"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}