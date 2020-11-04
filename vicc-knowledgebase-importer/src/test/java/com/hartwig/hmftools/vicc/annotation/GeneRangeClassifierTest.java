package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class GeneRangeClassifierTest {

    @Test
    public void canClassifyGeneLevelEvents() {
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("AKT1 act mut", "default_feature"));
    }

}