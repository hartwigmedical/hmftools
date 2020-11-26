package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class GeneLevelMatcherTest {

    private static final String EXON_KEYWORD = "exon";

    private static final Set<String> GENERIC_GENE_LEVEL_KEYWORDS = Sets.newHashSet("ALTERATION");
    private static final Set<String> INACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet("biallelic inactivation");
    private static final Set<String> ACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet("act mut", "positive");

    @Test
    public void canAssessWhetherEventIsGeneLevelEvent() {
        EventMatcher matcher = new GeneLevelMatcher(EXON_KEYWORD,
                GENERIC_GENE_LEVEL_KEYWORDS,
                ACTIVATING_GENE_LEVEL_KEYWORDS,
                INACTIVATING_GENE_LEVEL_KEYWORDS);

        assertTrue(matcher.matches("AKT1", "AKT1 act mut"));
        assertTrue(matcher.matches("AKT1", "AKT1"));
        assertTrue(matcher.matches("ATK1", "ALTERATION"));
        assertTrue(matcher.matches("TP53", "biallelic inactivation"));

        assertFalse(matcher.matches("NPM", "Exon 12 mutation"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}