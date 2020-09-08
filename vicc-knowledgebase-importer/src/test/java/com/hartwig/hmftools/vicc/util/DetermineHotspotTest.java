package com.hartwig.hmftools.vicc.util;

import static org.junit.Assert.*;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class DetermineHotspotTest {

    @Test
    public void canAssessWhetherFeatureIsProteinAnnotation() {
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("K5N"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("L2230V"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("V5del"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("L755_T759del"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("L755_T756delinsPP"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("D770delinsGY"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("G10dup"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("G10fs"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("G10fs*"));
        assertTrue(DetermineHotspot.isResolvableProteinAnnotation("*10L"));

        // Just plain wrong protein annotations
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation(Strings.EMPTY));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("truncating"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("20LtoV"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("L20"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("LP"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("L"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("L2"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("L20Pdel5"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("fs"));

        // Splice variants are ignored by hotspot extractor:
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("T599insTT"));

        // Inframe event is too long
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("L4_T40del"));
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("L698_S1037dup"));

        // Wild-type mutations are ignored
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(DetermineHotspot.isResolvableProteinAnnotation("V600E/K"));
    }

}