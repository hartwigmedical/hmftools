package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HotspotClassifierTest {

    @Test
    public void canAssessWhetherEventIsHotspot() {
        EventMatcher classifier = HotspotClassifier.create(Lists.newArrayList(), new ProteinAnnotationExtractor());

        assertTrue(classifier.matches("any", "V600E"));
    }

    @Test
    public void canAssessWhetherEventIsTypicalProteinAnnotation() {
        assertTrue(HotspotClassifier.isValidProteinAnnotation("K5N"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("L2230V"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("V5del"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("L755_T759del"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("L755_T756delinsPP"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("D770delinsGY"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("G10dup"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("G10fs"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("G10fs*"));
        assertTrue(HotspotClassifier.isValidProteinAnnotation("*10L"));

        // Just plain wrong protein annotations
        assertFalse(HotspotClassifier.isValidProteinAnnotation(Strings.EMPTY));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("truncating"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("20LtoV"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("L20"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("LP"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("L"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("L2"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("L20Pdel5"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("fs"));

        // Splice variants are ignored by hotspot extractor:
        assertFalse(HotspotClassifier.isValidProteinAnnotation("963_D1010splice"));

        // Frameshifts which also change the codon in which the frame is shifted are ignored.
        assertFalse(HotspotClassifier.isValidProteinAnnotation("G20Pfs"));

        // Frameshifts which are preceded by stop codon are ignored.
        assertFalse(HotspotClassifier.isValidProteinAnnotation("G20*fs"));

        // Not a correctly formatted insert (position not clear)
        assertFalse(HotspotClassifier.isValidProteinAnnotation("T599insTT"));

        // Not a correctly formatted insert
        assertFalse(HotspotClassifier.isValidProteinAnnotation("V600D_K601insFGLAT"));

        // Inframe event is too long
        assertFalse(HotspotClassifier.isValidProteinAnnotation("L4_T40del"));
        assertFalse(HotspotClassifier.isValidProteinAnnotation("L698_S1037dup"));

        // Wild-type mutations are ignored
        assertFalse(HotspotClassifier.isValidProteinAnnotation("V600X"));

        // Mutations with logical OR are ignored
        assertFalse(HotspotClassifier.isValidProteinAnnotation("V600E/K"));
    }
}