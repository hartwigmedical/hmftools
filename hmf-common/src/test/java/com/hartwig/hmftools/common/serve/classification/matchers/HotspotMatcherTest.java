package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HotspotMatcherTest {

    @Test
    public void canAssessWhetherEventIsHotspot() {
        EventMatcher matcher = testHotspotMatcher();

        assertTrue(matcher.matches("any", "K5N"));
        assertTrue(matcher.matches("any", "L2230V"));
        assertTrue(matcher.matches("any", "V5del"));
        assertTrue(matcher.matches("any", "L755_T759del"));
        assertTrue(matcher.matches("any", "L755_T756delinsPP"));
        assertTrue(matcher.matches("any", "D770delinsGY"));
        assertTrue(matcher.matches("any", "G10dup"));
        assertTrue(matcher.matches("any", "G10fs"));
        assertTrue(matcher.matches("any", "G10fs*"));
        assertTrue(matcher.matches("any", "*10L"));

        // Just plain wrong protein annotations
        assertFalse(matcher.matches("any", Strings.EMPTY));
        assertFalse(matcher.matches("any", "truncating"));
        assertFalse(matcher.matches("any", "20LtoV"));
        assertFalse(matcher.matches("any", "L20"));
        assertFalse(matcher.matches("any", "LP"));
        assertFalse(matcher.matches("any", "L"));
        assertFalse(matcher.matches("any", "L2"));
        assertFalse(matcher.matches("any", "L20Pdel5"));
        assertFalse(matcher.matches("any", "fs"));

        // Splice variants are ignored by hotspot classifier:
        assertFalse(matcher.matches("any", "963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(matcher.matches("any", "G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(matcher.matches("any", "G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(matcher.matches("any", "T599insTT"));

        // Not a correctly formatted insert
        assertFalse(matcher.matches("any", "V600D_K601insFGLAT"));

        // Wild-type mutations are ignored
        assertFalse(matcher.matches("any", "V600X"));

        // Mutations with logical OR are ignored
        assertFalse(matcher.matches("any", "V600E/K"));

        // Inframe event is too long
        assertFalse(matcher.matches("any", "L4_T40del"));
        assertFalse(matcher.matches("any", "L698_S1037dup"));

        // Ignore hotspots on fusion genes.
        assertFalse(matcher.matches("any", "EML4-ALK L1123R"));
    }

    @Test
    public void canDetermineComplexMatches() {
        HotspotMatcher matcher = testHotspotMatcher();

        assertFalse(matcher.isComplexMatch("any", "V600E"));
        assertFalse(matcher.isComplexMatch("any", "G10fs*"));
        assertFalse(matcher.isComplexMatch("any", "L755_T759del"));
        assertFalse(matcher.isComplexMatch("any", "EML4-ALK L1123R"));

        assertTrue(matcher.isComplexMatch("any", "L698_S1037dup"));
        assertTrue(matcher.isComplexMatch("any", "L4_T40del"));
    }

    @NotNull
    private static HotspotMatcher testHotspotMatcher() {
        return new HotspotMatcher(event -> event, new FusionPairMatcher(Sets.newHashSet(), Sets.newHashSet(), Sets.newHashSet()));
    }
}