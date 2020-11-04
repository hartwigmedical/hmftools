package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsHotspot() {
        assertTrue(HotspotClassifier.isHotspot("K5N"));
        assertTrue(HotspotClassifier.isHotspot("L2230V"));
        assertTrue(HotspotClassifier.isHotspot("V5del"));
        assertTrue(HotspotClassifier.isHotspot("L755_T759del"));
        assertTrue(HotspotClassifier.isHotspot("L755_T756delinsPP"));
        assertTrue(HotspotClassifier.isHotspot("D770delinsGY"));
        assertTrue(HotspotClassifier.isHotspot("G10dup"));
        assertTrue(HotspotClassifier.isHotspot("G10fs"));
        assertTrue(HotspotClassifier.isHotspot("G10fs*"));
        assertTrue(HotspotClassifier.isHotspot("*10L"));

        // Just plain wrong protein annotations
        assertFalse(HotspotClassifier.isHotspot(Strings.EMPTY));
        assertFalse(HotspotClassifier.isHotspot("truncating"));
        assertFalse(HotspotClassifier.isHotspot("20LtoV"));
        assertFalse(HotspotClassifier.isHotspot("L20"));
        assertFalse(HotspotClassifier.isHotspot("LP"));
        assertFalse(HotspotClassifier.isHotspot("L"));
        assertFalse(HotspotClassifier.isHotspot("L2"));
        assertFalse(HotspotClassifier.isHotspot("L20Pdel5"));
        assertFalse(HotspotClassifier.isHotspot("fs"));

        // Hotspots on fusion genes are ignored
        assertFalse(HotspotClassifier.isHotspot("EML4-ALK L1152R"));

        // Splice variants are ignored by hotspot extractor:
        assertFalse(HotspotClassifier.isHotspot("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(HotspotClassifier.isHotspot("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(HotspotClassifier.isHotspot("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(HotspotClassifier.isHotspot("T599insTT"));

        // Not a correctly formatted insert
        assertFalse(HotspotClassifier.isHotspot("V600D_K601insFGLAT"));

        // Inframe event is too long
        assertFalse(HotspotClassifier.isHotspot("L4_T40del"));
        assertFalse(HotspotClassifier.isHotspot("L698_S1037dup"));

        // Wild-type mutations are ignored
        assertFalse(HotspotClassifier.isHotspot("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(HotspotClassifier.isHotspot("V600E/K"));
    }

    @Test
    public void canConvertFeatureNameToProteinAnnotation() {
        assertEquals("E709K", HotspotClassifier.extractProteinAnnotation("E709K"));
        assertEquals("E709K", HotspotClassifier.extractProteinAnnotation("EGFR E709K "));
        assertEquals("E709K", HotspotClassifier.extractProteinAnnotation("EGFR:E709K"));
        assertEquals("E709K", HotspotClassifier.extractProteinAnnotation("EGFR:p.E709K"));
        assertEquals("E709K", HotspotClassifier.extractProteinAnnotation("EGFR p.E709K"));
        assertEquals("E709K", HotspotClassifier.extractProteinAnnotation("E709K (c.2100A>c)"));

        assertEquals("G778_P780dup", HotspotClassifier.extractProteinAnnotation("G778_P780DUP"));
        assertEquals("V560del", HotspotClassifier.extractProteinAnnotation("KIT:p.V560DEL"));
        assertEquals("V560fs", HotspotClassifier.extractProteinAnnotation("KIT:p.V560FS"));
        assertEquals("V560fs", HotspotClassifier.extractProteinAnnotation("KIT:p.V560FS*"));
        assertEquals("V560insAYVM", HotspotClassifier.extractProteinAnnotation("KIT:p.V560INSAYVM"));
        assertEquals("V560insINS", HotspotClassifier.extractProteinAnnotation("KIT:p.V560INSINS"));
        assertEquals("V560delinsDEL", HotspotClassifier.extractProteinAnnotation("KIT:p.V560DELINSDEL"));
    }
}