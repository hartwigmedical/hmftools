package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CombinedClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsFusionAndExonRange() {
        assertTrue(CombinedClassifier.isFusionPairAndGeneRangeExon("EXON 14 SKIPPING MUTATION", "MET"));
        assertFalse(CombinedClassifier.isFusionPairAndGeneRangeExon("EXON 14 SKIPPING MUTATION", "NRG1"));
    }

    @Test
    public void canAssessWhetherFeatureIsCombinedEvent() {
        assertTrue(CombinedClassifier.isCombinedEvent("Ex19 del L858R", "EGFR"));

        assertTrue(CombinedClassifier.isCombinedEvent("AR (F877L) + AR (T878A)", "AR"));
        assertTrue(CombinedClassifier.isCombinedEvent("JAK1 (S646F;R683)", "JAK1"));
        assertTrue(CombinedClassifier.isCombinedEvent("KIT (627-664,664-714,449-514)", "KIT"));
        assertFalse(CombinedClassifier.isCombinedEvent("KIT mutation in exon 9,11,13,14 or 17", "KIT"));

        assertTrue(CombinedClassifier.isCombinedEvent("BCR-ABL F486S", "ABL"));
        assertTrue(CombinedClassifier.isCombinedEvent("NPM1-ALK  amp", "ALK"));

        assertFalse(CombinedClassifier.isCombinedEvent("Exon 20 insertions/deletions", "ERBB2"));
        assertFalse(CombinedClassifier.isCombinedEvent("3'UTR alteration (c.642+70C>A)", "VHL"));
    }
}