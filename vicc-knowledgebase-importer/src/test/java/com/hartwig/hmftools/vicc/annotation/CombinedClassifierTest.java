package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CombinedClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsCombinedEvent() {
        assertTrue(CombinedClassifier.isCombinedEvent("Ex19 del L858R", "EGFR"));
        assertTrue(CombinedClassifier.isCombinedEvent("BCR-ABL F486S", "ABL"));
        assertTrue(CombinedClassifier.isCombinedEvent("NPM1-ALK  amp", "ALK"));
    }
}