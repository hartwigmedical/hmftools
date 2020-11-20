package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.junit.Test;

public class AmplificationClassifierTest {

    @Test
    public void canAssessWhetherEventIsAmplification() {
        EventMatcher classifier = AmplificationClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("ALK", "ALK over exp"));
        assertTrue(classifier.matches("ALK", "ALK  amp"));

        assertFalse(classifier.matches("BRAF", "V600E"));
    }
}