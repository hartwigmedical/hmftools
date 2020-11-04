package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CombinedClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsCombinedEvent() {
        assertTrue(CombinedClassifier.isCombinedEvent("BCR-ABL F486S"));
        assertTrue(CombinedClassifier.isCombinedEvent("NPM1-ALK amp"));
    }
}