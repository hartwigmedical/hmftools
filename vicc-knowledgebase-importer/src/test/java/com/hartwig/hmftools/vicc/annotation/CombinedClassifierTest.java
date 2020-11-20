package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class CombinedClassifierTest {

    @Test
    public void canAssessWhetherEventIsCombinedEvent() {
        EventMatcher classifier = new CombinedClassifier();

        assertTrue(classifier.matches("EGFR", "Ex19 del L858R"));

        assertTrue(classifier.matches("AR", "AR (F877L) + AR (T878A)"));
        assertTrue(classifier.matches("JAK1", "JAK1 (S646F;R683)"));
        assertTrue(classifier.matches("KIT", "KIT (627-664,664-714,449-514)"));
        assertFalse(classifier.matches("KIT", "KIT mutation in exon 9,11,13,14 or 17"));

        assertTrue(classifier.matches("ABL", "BCR-ABL F486S"));
        assertTrue(classifier.matches("ALK", "NPM1-ALK  amp"));

        assertFalse(classifier.matches("ERBB2", "Exon 20 insertions/deletions"));
        assertFalse(classifier.matches("VHL", "3'UTR alteration (c.642+70C>A)"));
    }
}