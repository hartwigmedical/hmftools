package com.hartwig.hmftools.common.serve.classification;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotClassifierTest {

    @Test
    public void canAssessWhetherEventIsHotspot() {
        EventMatcher classifier = HotspotClassifier.create(Lists.newArrayList(), event -> event);

        assertTrue(classifier.matches("any", "V600E"));
        assertTrue(classifier.matches("any", "V5del"));

        // Inframe event is too long -> still protein annotation.
        assertFalse(classifier.matches("any", "L4_T40del"));
        assertFalse(classifier.matches("any", "L698_S1037dup"));

        // Do not consider hotspots on fusion genes.
        assertFalse(classifier.matches("ALK", "EML4-ALK L100R"));
    }

    @Test
    public void canAssessWhetherEventIsProteinAnnotation() {
        assertTrue(HotspotClassifier.isProteinAnnotation("K5N"));
        assertTrue(HotspotClassifier.isProteinAnnotation("L2230V"));
        assertTrue(HotspotClassifier.isProteinAnnotation("V5del"));
        assertTrue(HotspotClassifier.isProteinAnnotation("L755_T759del"));
        assertTrue(HotspotClassifier.isProteinAnnotation("L4_T40del"));
        assertTrue(HotspotClassifier.isProteinAnnotation("L698_S1037dup"));
        assertTrue(HotspotClassifier.isProteinAnnotation("L755_T756delinsPP"));
        assertTrue(HotspotClassifier.isProteinAnnotation("D770delinsGY"));
        assertTrue(HotspotClassifier.isProteinAnnotation("G10dup"));
        assertTrue(HotspotClassifier.isProteinAnnotation("G10fs"));
        assertTrue(HotspotClassifier.isProteinAnnotation("G10fs*"));
        assertTrue(HotspotClassifier.isProteinAnnotation("*10L"));

        // Just plain wrong protein annotations
        assertFalse(HotspotClassifier.isProteinAnnotation(Strings.EMPTY));
        assertFalse(HotspotClassifier.isProteinAnnotation("truncating"));
        assertFalse(HotspotClassifier.isProteinAnnotation("20LtoV"));
        assertFalse(HotspotClassifier.isProteinAnnotation("L20"));
        assertFalse(HotspotClassifier.isProteinAnnotation("LP"));
        assertFalse(HotspotClassifier.isProteinAnnotation("L"));
        assertFalse(HotspotClassifier.isProteinAnnotation("L2"));
        assertFalse(HotspotClassifier.isProteinAnnotation("L20Pdel5"));
        assertFalse(HotspotClassifier.isProteinAnnotation("fs"));

        // Splice variants are ignored by hotspot classifier:
        assertFalse(HotspotClassifier.isProteinAnnotation("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(HotspotClassifier.isProteinAnnotation("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(HotspotClassifier.isProteinAnnotation("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(HotspotClassifier.isProteinAnnotation("T599insTT"));

        // Not a correctly formatted insert
        assertFalse(HotspotClassifier.isProteinAnnotation("V600D_K601insFGLAT"));

        // Wild-type mutations are ignored
        assertFalse(HotspotClassifier.isProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(HotspotClassifier.isProteinAnnotation("V600E/K"));
    }
}