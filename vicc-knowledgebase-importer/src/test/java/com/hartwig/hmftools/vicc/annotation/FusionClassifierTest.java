package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FusionClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsFusionPair() {
        assertTrue(FusionClassifier.isFusionPair("ABL1-BCR fusion", "ABL1"));
    }

    @Test
    public void canAssessWhetherFeatureIsPromiscuousFusion() {
        assertTrue(FusionClassifier.isPromiscuousFusion("BRAF fusion", "BRAF"));
    }

}