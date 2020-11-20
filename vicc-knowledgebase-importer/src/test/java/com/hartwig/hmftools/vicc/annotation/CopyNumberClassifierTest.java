package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CopyNumberClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsAmplification() {
        assertTrue(CopyNumberClassifier.isAmplification("ALK", "ALK over exp"));
        assertTrue(CopyNumberClassifier.isAmplification("ALK", "ALK  amp"));

        assertFalse(CopyNumberClassifier.isAmplification("MET", "MET amplification + mutation"));
        assertFalse(CopyNumberClassifier.isAmplification("ALK", "NPM1-ALK  amp"));
        assertFalse(CopyNumberClassifier.isAmplification("BRAF", "V600E"));
    }

    @Test
    public void canAssessWhetherFeatureIsDeletion() {
        assertTrue(CopyNumberClassifier.isDeletion("CDKN2A", "CDKN2A del"));
        assertTrue(CopyNumberClassifier.isDeletion("CDKN2A", "CDKN2A dec exp"));

        assertFalse(CopyNumberClassifier.isDeletion("MET", "MET deletion + mutation"));
        assertFalse(CopyNumberClassifier.isDeletion("EGFR", "EGFR inframe deletion"));
        assertFalse(CopyNumberClassifier.isDeletion("BRAF", "V600E"));
    }
}