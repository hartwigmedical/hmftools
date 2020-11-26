package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.junit.Test;

public class FusionPairAndExonRangeMatcherTest {

    private static final Map<String, Set<String>> FUSION_PAIR_AND_EXON_RANGES_PER_GENE = Maps.newHashMap();

    static {
        FUSION_PAIR_AND_EXON_RANGES_PER_GENE.put("MET", Sets.newHashSet("EXON 14 SKIPPING MUTATION"));
    }

    @Test
    public void canAssessWhetherEventIsFusionAndExonRange() {
        EventMatcher matcher = new FusionPairAndExonRangeMatcher(FUSION_PAIR_AND_EXON_RANGES_PER_GENE);

        assertTrue(matcher.matches("MET", "EXON 14 SKIPPING MUTATION"));
        assertFalse(matcher.matches("NRG1", "EXON 14 SKIPPING MUTATION"));
    }
}