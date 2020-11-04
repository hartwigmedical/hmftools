package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CopyNumberClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsCopyNumber() {
        assertTrue(CopyNumberClassifier.isAmplification("ALK  amp", "xx"));
        assertTrue(CopyNumberClassifier.isDeletion("CDKN2A del", "xx"));

        assertFalse(CopyNumberClassifier.isAmplification("NPM1-ALK  amp", "xx"));
    }
}