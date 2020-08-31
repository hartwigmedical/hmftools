package com.hartwig.hmftools.serve.hotspot;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotGeneratorTest {

    @Test
    public void canAssessWhetherFeatureIsProteinAnnotation() {
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("K5N"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("L2230V"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("V5del"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("L755_T759del"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("L755_T756delinsPP"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("D770delinsGY"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("G10dup"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("G10fs"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("G10fs*"));
        assertTrue(HotspotGenerator.isResolvableProteinAnnotation("*10L"));

        // Just plain wrong protein annotations
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation(Strings.EMPTY));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("truncating"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("20LtoV"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("L20"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("LP"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("L"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("L2"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("L20Pdel5"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("fs"));

        // Splice variants are ignored by hotspot extractor:
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("T599insTT"));

        // Inframe event is too long
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("L4_T40del"));
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("L698_S1037dup"));

        // Wild-type mutations are ignored
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(HotspotGenerator.isResolvableProteinAnnotation("V600E/K"));
    }
}