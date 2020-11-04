package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FusionClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsFusionPair() {
        assertTrue(FusionClassifier.isFusionPair("ABL1-BCR fusion", "ABL1"));
        assertTrue(FusionClassifier.isFusionPair("EML4-ALK", "ALK"));
        assertTrue(FusionClassifier.isFusionPair("VIII", "EGFR"));

        assertFalse(FusionClassifier.isFusionPair("p61BRAF-V600E", "BRAF"));
        assertFalse(FusionClassifier.isFusionPair("BRAF fusion + mutation", "BRAF"));
    }

    @Test
    public void canAssessWhetherFeatureIsPromiscuousFusion() {
        assertTrue(FusionClassifier.isPromiscuousFusion("BRAF fusion", "BRAF"));
        assertTrue(FusionClassifier.isPromiscuousFusion("ROS1 REARRANGEMENT", "ROS1"));

        assertFalse(FusionClassifier.isPromiscuousFusion("BRAF fusion + mutation", "BRAF"));
    }

}