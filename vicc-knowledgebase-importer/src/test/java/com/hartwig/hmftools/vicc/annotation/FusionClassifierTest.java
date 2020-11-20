package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FusionClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsFusionPair() {
        assertTrue(FusionClassifier.isFusionPair("ABL1", "ABL1-BCR fusion"));
        assertTrue(FusionClassifier.isFusionPair("ALK", "EML4-ALK"));
        assertTrue(FusionClassifier.isFusionPair("EGFR", "VIII"));
        assertTrue(FusionClassifier.isFusionPair("NKX2-1", "TRB-NKX2-1 Fusion"));

        assertFalse(FusionClassifier.isFusionPair("BRAF", "BRAF amp"));
        assertFalse(FusionClassifier.isFusionPair("AR", "AR-V7"));
        assertFalse(FusionClassifier.isFusionPair("BRAF", "BRAF fusion + mutation"));
        assertFalse(FusionClassifier.isFusionPair("BRAF", "V600E"));
    }

    @Test
    public void canAssessWhetherFeatureIsPromiscuousFusion() {
        assertTrue(FusionClassifier.isPromiscuousFusion("BRAF", "BRAF fusion"));
        assertTrue(FusionClassifier.isPromiscuousFusion("ROS1", "ROS1 REARRANGEMENT"));

        assertFalse(FusionClassifier.isPromiscuousFusion("BRAF", "BRAF fusion + mutation"));
        assertFalse(FusionClassifier.isPromiscuousFusion("BRAF", "V600E"));
    }
}