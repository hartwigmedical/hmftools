package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class SignatureClassifierTest {

    @Test
    public void canAssessWhetherEventIsSignature() {
        EventMatcher classifier = SignatureClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("-", "Microsatellite Instability-High"));

        assertFalse(classifier.matches("BRAF", "V600E"));
    }
}