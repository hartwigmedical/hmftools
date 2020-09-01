package com.hartwig.hmftools.serve.hotspot;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ProteinToHotspotConverterTest {

    @Test
    public void canAssessWhetherFeatureIsProteinAnnotation() {
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("K5N"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("L2230V"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("V5del"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("L755_T759del"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("L755_T756delinsPP"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("D770delinsGY"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("G10dup"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("G10fs"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("G10fs*"));
        assertTrue(ProteinToHotspotConverter.isResolvableProteinAnnotation("*10L"));

        // Just plain wrong protein annotations
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation(Strings.EMPTY));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("truncating"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("20LtoV"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("L20"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("LP"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("L"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("L2"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("L20Pdel5"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("fs"));

        // Splice variants are ignored by hotspot extractor:
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("T599insTT"));

        // Inframe event is too long
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("L4_T40del"));
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("L698_S1037dup"));

        // Wild-type mutations are ignored
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(ProteinToHotspotConverter.isResolvableProteinAnnotation("V600E/K"));
    }
}