package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.junit.Test;

public class FusionPairAndExonMatcherTest {

    private static final Map<String, Set<String>> FUSION_PAIR_AND_EXONS_PER_GENE = Maps.newHashMap();

    static {
        FUSION_PAIR_AND_EXONS_PER_GENE.put("MET", Sets.newHashSet("EXON 14 SKIPPING MUTATION"));
    }

    @Test
    public void canAssessWhetherEventIsFusionAndExon() {
        EventMatcher matcher = new FusionPairAndExonMatcher(FUSION_PAIR_AND_EXONS_PER_GENE);

        assertTrue(matcher.matches("MET", "EXON 14 SKIPPING MUTATION"));
        assertFalse(matcher.matches("NRG1", "EXON 14 SKIPPING MUTATION"));
    }
}