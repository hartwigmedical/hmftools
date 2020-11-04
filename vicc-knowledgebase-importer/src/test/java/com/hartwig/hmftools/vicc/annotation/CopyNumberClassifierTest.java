package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CopyNumberClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsAmplification() {
        assertTrue(CopyNumberClassifier.isAmplification("ALK over exp", "ALK"));
        assertTrue(CopyNumberClassifier.isAmplification("ALK  amp", "ALK"));

        assertFalse(CopyNumberClassifier.isAmplification("MET amplification + mutation", "MET"));
        assertFalse(CopyNumberClassifier.isAmplification("NPM1-ALK  amp", "ALK"));
    }

    @Test
    public void canAssessWhetherFeatureIsDeletion() {
        assertTrue(CopyNumberClassifier.isDeletion("CDKN2A del", "CDKN2A"));

        assertFalse(CopyNumberClassifier.isDeletion("MET deletion + mutation", "MET"));
        assertFalse(CopyNumberClassifier.isDeletion("EGFR inframe deletion", "EGFR"));
    }
}