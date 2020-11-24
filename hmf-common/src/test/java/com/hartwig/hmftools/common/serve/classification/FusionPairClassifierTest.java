package com.hartwig.hmftools.common.serve.classification;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.junit.Test;

public class FusionPairClassifierTest {

    @Test
    public void canAssessWhetherEventIsFusionPair() {
        EventMatcher classifier = FusionPairClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("ABL1", "ABL1-BCR fusion"));
        assertTrue(classifier.matches("ALK", "EML4-ALK"));
        assertTrue(classifier.matches("EGFR", "VIII"));
        assertTrue(classifier.matches("NKX2-1", "TRB-NKX2-1 Fusion"));

        assertFalse(classifier.matches("BRAF", "BRAF amp"));
        assertFalse(classifier.matches("AR", "AR-V7"));
        assertFalse(classifier.matches("BRAF", "V600E"));
    }
}