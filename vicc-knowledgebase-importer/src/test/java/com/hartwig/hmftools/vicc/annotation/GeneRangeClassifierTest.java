package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class GeneRangeClassifierTest {

    @Test
    public void canClassifyGeneLevelEvents() {
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("AKT1 act mut", "default_feature"));
    }

    @Test
    public void canClassifyGeneRangeCodonEvents() {
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("GNAS (R201)"));
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("E709X"));
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("EGFR E709X "));
    }
}