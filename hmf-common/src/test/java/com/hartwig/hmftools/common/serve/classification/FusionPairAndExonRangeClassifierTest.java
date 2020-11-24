package com.hartwig.hmftools.common.serve.classification;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FusionPairAndExonRangeClassifierTest {

    @Test
    public void canAssessWhetherEventIsFusionAndExonRange() {
        EventMatcher classifier = new FusionPairAndExonRangeClassifier();

        assertTrue(classifier.matches("MET", "EXON 14 SKIPPING MUTATION"));
        assertFalse(classifier.matches("NRG1", "EXON 14 SKIPPING MUTATION"));
    }
}