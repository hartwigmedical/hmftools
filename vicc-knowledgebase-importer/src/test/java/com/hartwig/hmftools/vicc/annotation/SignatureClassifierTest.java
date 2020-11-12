package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SignatureClassifierTest {

    @Test
    public void canAssessWhetherFeatureIsSignature() {
        assertTrue(SignatureClassifier.isSignature("Microsatellite Instability-High"));

        assertFalse(SignatureClassifier.isSignature("V600E"));
    }
}