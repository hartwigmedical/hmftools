package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FusionClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsFusion() {
        assertTrue(FusionClassifier.isFusionPair("ABL1-BCR fusion", "ABL1", "splice"));
    }

}