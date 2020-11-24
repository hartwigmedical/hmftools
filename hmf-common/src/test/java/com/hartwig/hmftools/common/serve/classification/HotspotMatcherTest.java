package com.hartwig.hmftools.common.serve.classification;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotMatcherTest {

    @Test
    public void canAssessWhetherEventIsHotspot() {
        EventMatcher matcher = new HotspotMatcher(event -> event);

        assertTrue(matcher.matches("any", "V600E"));
        assertTrue(matcher.matches("any", "V5del"));

        // Inframe event is too long -> still protein annotation.
        assertFalse(matcher.matches("any", "L4_T40del"));
        assertFalse(matcher.matches("any", "L698_S1037dup"));

        // Do not consider hotspots on fusion genes.
        assertFalse(matcher.matches("ALK", "EML4-ALK L100R"));
    }

    @Test
    public void canAssessWhetherEventIsProteinAnnotation() {
        assertTrue(HotspotMatcher.isProteinAnnotation("K5N"));
        assertTrue(HotspotMatcher.isProteinAnnotation("L2230V"));
        assertTrue(HotspotMatcher.isProteinAnnotation("V5del"));
        assertTrue(HotspotMatcher.isProteinAnnotation("L755_T759del"));
        assertTrue(HotspotMatcher.isProteinAnnotation("L4_T40del"));
        assertTrue(HotspotMatcher.isProteinAnnotation("L698_S1037dup"));
        assertTrue(HotspotMatcher.isProteinAnnotation("L755_T756delinsPP"));
        assertTrue(HotspotMatcher.isProteinAnnotation("D770delinsGY"));
        assertTrue(HotspotMatcher.isProteinAnnotation("G10dup"));
        assertTrue(HotspotMatcher.isProteinAnnotation("G10fs"));
        assertTrue(HotspotMatcher.isProteinAnnotation("G10fs*"));
        assertTrue(HotspotMatcher.isProteinAnnotation("*10L"));

        // Just plain wrong protein annotations
        assertFalse(HotspotMatcher.isProteinAnnotation(Strings.EMPTY));
        assertFalse(HotspotMatcher.isProteinAnnotation("truncating"));
        assertFalse(HotspotMatcher.isProteinAnnotation("20LtoV"));
        assertFalse(HotspotMatcher.isProteinAnnotation("L20"));
        assertFalse(HotspotMatcher.isProteinAnnotation("LP"));
        assertFalse(HotspotMatcher.isProteinAnnotation("L"));
        assertFalse(HotspotMatcher.isProteinAnnotation("L2"));
        assertFalse(HotspotMatcher.isProteinAnnotation("L20Pdel5"));
        assertFalse(HotspotMatcher.isProteinAnnotation("fs"));

        // Splice variants are ignored by hotspot classifier:
        assertFalse(HotspotMatcher.isProteinAnnotation("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(HotspotMatcher.isProteinAnnotation("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(HotspotMatcher.isProteinAnnotation("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(HotspotMatcher.isProteinAnnotation("T599insTT"));

        // Not a correctly formatted insert
        assertFalse(HotspotMatcher.isProteinAnnotation("V600D_K601insFGLAT"));

        // Wild-type mutations are ignored
        assertFalse(HotspotMatcher.isProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(HotspotMatcher.isProteinAnnotation("V600E/K"));
    }
}