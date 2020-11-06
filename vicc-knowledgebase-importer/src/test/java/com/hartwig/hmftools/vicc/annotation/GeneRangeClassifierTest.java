package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class GeneRangeClassifierTest {

    @Test
    public void canClassifyGeneLevelEvents() {
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("AKT1 act mut", "AKT1"));
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("AKT1", "AKT1"));
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("biallelic inactivation", "TP53"));
    }

    @Test
    public void canClassifyGeneRangeCodonEvents() {
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("GNAS (R201)"));
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("E709X"));
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("EGFR E709X "));

        assertFalse(GeneRangeClassifier.isGeneRangeCodonEvent("EZH2 (Y641,A677)"));
    }
}