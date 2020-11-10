package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class GeneRangeClassifierTest {

    @Test
    public void canClassifyGeneLevelEvents() {
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("AKT1 act mut", "AKT1"));
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("AKT1", "AKT1"));
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("positive", "AKT1"));
        assertTrue(GeneRangeClassifier.isGeneLevelEvent("biallelic inactivation", "TP53"));

        assertFalse(GeneRangeClassifier.isGeneLevelEvent("3' UTR MUTATION", "KIT"));
        assertFalse(GeneRangeClassifier.isGeneLevelEvent("Exon 12 mutation", "NPM"));
        assertFalse(GeneRangeClassifier.isGeneLevelEvent("V600E", "BRAF"));
    }

    @Test
    public void canClassifyGeneRangeExonEvents() {
        assertTrue(GeneRangeClassifier.isGeneRangeExonEvent("EGFR exon 19 deletions", "EGFR"));
        assertTrue(GeneRangeClassifier.isGeneRangeExonEvent("EXON 12 MUTATION", "AXSL1"));

        assertFalse(GeneRangeClassifier.isGeneRangeExonEvent("EGFR exon 19 deletions + amp", "EGFR"));
        assertFalse(GeneRangeClassifier.isGeneRangeExonEvent("V600E", "BRAF"));
    }

    @Test
    public void canClassifyGeneRangeCodonEvents() {
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("GNAS (R201)"));
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("E709X"));
        assertTrue(GeneRangeClassifier.isGeneRangeCodonEvent("EGFR E709X "));

        assertFalse(GeneRangeClassifier.isGeneRangeCodonEvent(Strings.EMPTY));
        assertFalse(GeneRangeClassifier.isGeneRangeCodonEvent("600"));
        assertFalse(GeneRangeClassifier.isGeneRangeCodonEvent("EZH2 (Y641,A677)"));
        assertFalse(GeneRangeClassifier.isGeneRangeCodonEvent("T148HFSX9"));
        assertFalse(GeneRangeClassifier.isGeneRangeCodonEvent("RARE EX 18-21 MUT"));
        assertFalse(GeneRangeClassifier.isGeneRangeCodonEvent("V600E"));
    }

    @Test
    public void canCountDigitSequences() {
        assertEquals(0, GeneRangeClassifier.countDigitSequences(Strings.EMPTY));
        assertEquals(0, GeneRangeClassifier.countDigitSequences("hi"));
        assertEquals(1, GeneRangeClassifier.countDigitSequences("V600K"));

        assertEquals(4, GeneRangeClassifier.countDigitSequences("A1B2C3D4"));
        assertEquals(2, GeneRangeClassifier.countDigitSequences("100A200"));
    }
}