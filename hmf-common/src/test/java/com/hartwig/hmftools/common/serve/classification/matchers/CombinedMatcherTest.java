package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CombinedMatcherTest {

    @Test
    public void canAssessWhetherEventIsCombinedEvent() {
        EventMatcher matcher = buildCombinedMatcher();

        assertTrue(matcher.matches("EGFR", "Ex19 del L858R"));

        assertTrue(matcher.matches("AR", "AR (F877L) + AR (T878A)"));
        assertTrue(matcher.matches("JAK1", "JAK1 (S646F;R683)"));
        assertTrue(matcher.matches("KIT", "KIT (627-664,664-714,449-514)"));
        assertFalse(matcher.matches("KIT", "KIT mutation in exon 9,11,13,14 or 17"));

        assertTrue(matcher.matches("ABL", "BCR-ABL F486S"));
        assertTrue(matcher.matches("ALK", "NPM1-ALK  amp"));

        assertFalse(matcher.matches("ERBB2", "Exon 20 insertions/deletions"));
        assertFalse(matcher.matches("VHL", "3'UTR alteration (c.642+70C>A)"));
    }

    @NotNull
    private static EventMatcher buildCombinedMatcher() {
        FusionPairMatcher fusionPairMatcher = new FusionPairMatcher();
        HotspotMatcher hotspotMatcher = new HotspotMatcher(event -> event, fusionPairMatcher);

        return new CombinedMatcher(hotspotMatcher, fusionPairMatcher, new AmplificationMatcher());
    }
}