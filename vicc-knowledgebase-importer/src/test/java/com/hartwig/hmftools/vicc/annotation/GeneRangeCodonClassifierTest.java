package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class GeneRangeCodonClassifierTest {

    @Test
    public void canAssessWhetherEventIsGeneRangeCodon() {
        EventMatcher classifier = GeneRangeCodonClassifier.create(Lists.newArrayList());

        assertTrue(classifier.matches("GNAS", "GNAS (R201)"));
        assertTrue(classifier.matches("EGFR", "E709X"));
        assertTrue(classifier.matches("EGFR", "EGFR E709X "));

        assertFalse(classifier.matches(Strings.EMPTY, Strings.EMPTY));
        assertFalse(classifier.matches("BRAF", "600"));
        assertFalse(classifier.matches("EZH2", "EZH2 (Y641,A677)"));
        assertFalse(classifier.matches("EZH2", "T148HFSX9"));
        assertFalse(classifier.matches("EGFR", "RARE EX 18-21 MUT"));
        assertFalse(classifier.matches("BRAF", "V600E"));
    }

    @Test
    public void canCountDigitSequences() {
        assertEquals(0, GeneRangeCodonClassifier.countDigitSequences(Strings.EMPTY));
        assertEquals(0, GeneRangeCodonClassifier.countDigitSequences("hi"));
        assertEquals(1, GeneRangeCodonClassifier.countDigitSequences("V600K"));

        assertEquals(4, GeneRangeCodonClassifier.countDigitSequences("A1B2C3D4"));
        assertEquals(2, GeneRangeCodonClassifier.countDigitSequences("100A200"));
    }
}