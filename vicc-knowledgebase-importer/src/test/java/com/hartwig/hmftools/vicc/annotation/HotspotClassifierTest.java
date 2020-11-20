package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotClassifierTest {

    @Test
    public void canAssessWhetherEventIsHotspot() {
        EventClassifier classifier = HotspotClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("any", "V600E"));
    }

    @Test
    public void canAssessWhetherEventIsTypicalProteinAnnotation() {
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("K5N"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("L2230V"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("V5del"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("L755_T759del"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("L755_T756delinsPP"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("D770delinsGY"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("G10dup"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("G10fs"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("G10fs*"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("G10fs*4"));
        assertTrue(HotspotClassifier.isTypicalProteinAnnotation("*10L"));

        // Just plain wrong protein annotations
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation(Strings.EMPTY));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("truncating"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("20LtoV"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("L20"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("LP"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("L"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("L2"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("L20Pdel5"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("fs"));

        // Splice variants are ignored by hotspot extractor:
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("T599insTT"));

        // Not a correctly formatted insert
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("V600D_K601insFGLAT"));

        // Inframe event is too long
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("L4_T40del"));
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("L698_S1037dup"));

        // Wild-type mutations are ignored
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(HotspotClassifier.isTypicalProteinAnnotation("V600E/K"));
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

        assertEquals("S978fs", HotspotClassifier.extractProteinAnnotation("S978FS*4"));
        assertEquals("S978fs", HotspotClassifier.extractProteinAnnotation("S978fs*123"));
    }
}