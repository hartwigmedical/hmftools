package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CombinedClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsFusionAndExonRange() {
        assertTrue(CombinedClassifier.isFusionPairAndGeneRangeExon("MET", "EXON 14 SKIPPING MUTATION"));
        assertFalse(CombinedClassifier.isFusionPairAndGeneRangeExon("NRG1", "EXON 14 SKIPPING MUTATION"));
    }

    @Test
    public void canAssessWhetherFeatureIsCombinedEvent() {
        assertTrue(CombinedClassifier.isCombinedEvent("EGFR", "Ex19 del L858R"));

        assertTrue(CombinedClassifier.isCombinedEvent("AR", "AR (F877L) + AR (T878A)"));
        assertTrue(CombinedClassifier.isCombinedEvent("JAK1", "JAK1 (S646F;R683)"));
        assertTrue(CombinedClassifier.isCombinedEvent("KIT", "KIT (627-664,664-714,449-514)"));
        assertFalse(CombinedClassifier.isCombinedEvent("KIT", "KIT mutation in exon 9,11,13,14 or 17"));

        assertTrue(CombinedClassifier.isCombinedEvent("ABL", "BCR-ABL F486S"));
        assertTrue(CombinedClassifier.isCombinedEvent("ALK", "NPM1-ALK  amp"));

        assertFalse(CombinedClassifier.isCombinedEvent("ERBB2", "Exon 20 insertions/deletions"));
        assertFalse(CombinedClassifier.isCombinedEvent("VHL", "3'UTR alteration (c.642+70C>A)"));
    }
}