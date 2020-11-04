package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CopyNumberClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsAmplification() {
        assertTrue(CopyNumberClassifier.isAmplification("ALK over exp"));
        assertTrue(CopyNumberClassifier.isAmplification("ALK  amp"));

        assertFalse(CopyNumberClassifier.isAmplification("MET amplification + mutation"));
        assertFalse(CopyNumberClassifier.isAmplification("NPM1-ALK  amp"));
    }

    @Test
    public void canAssessWhetherFeatureIsDeletion() {
        assertTrue(CopyNumberClassifier.isDeletion("CDKN2A del"));

        assertFalse(CopyNumberClassifier.isDeletion("MET deletion + mutation"));
        assertFalse(CopyNumberClassifier.isDeletion("EGFR inframe deletion"));
    }
}