package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class GeneLevelMatcherTest {

    private static final Set<String> BLACKLIST_KEY_PHRASES = Sets.newHashSet("Exon");
    private static final Set<String> GENERIC_GENE_LEVEL_KEY_PHRASES = Sets.newHashSet("ALTERATION");
    private static final Set<String> INACTIVATING_GENE_LEVEL_KEY_PHRASES = Sets.newHashSet("biallelic inactivation", "inact mut");
    private static final Set<String> ACTIVATING_GENE_LEVEL_KEY_PHRASES = Sets.newHashSet("act mut");

    @Test
    public void canAssessWhetherEventIsGeneLevelEvent() {
        EventMatcher matcher = new GeneLevelMatcher(BLACKLIST_KEY_PHRASES,
                GENERIC_GENE_LEVEL_KEY_PHRASES,
                ACTIVATING_GENE_LEVEL_KEY_PHRASES,
                INACTIVATING_GENE_LEVEL_KEY_PHRASES);

        assertTrue(matcher.matches("AKT1", "AKT1 act mut"));
        assertTrue(matcher.matches("AKT1", "AKT1"));
        assertTrue(matcher.matches("ATK1", "ALTERATION"));
        assertTrue(matcher.matches("TP53", "biallelic inactivation"));

        assertFalse(matcher.matches("NPM", "Exon 12 ALTERATION"));
        assertFalse(matcher.matches("BRAF", "V600E"));
    }
}