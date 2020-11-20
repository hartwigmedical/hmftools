package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ComplexClassifierTest {

    @Test
    public void canAssessWhetherEventIsComplexEvent() {
        assertTrue(ComplexClassifier.isComplexEvent("VHL", "Splicing alteration (c.464-2A>G)"));

        assertTrue(ComplexClassifier.isComplexEvent("KRAS", "KRAS ."));
        assertTrue(ComplexClassifier.isComplexEvent("APC", "APC p.I1557*fs*1"));
        assertTrue(ComplexClassifier.isComplexEvent("BRCA1", "BRCA1 L631Qfs*4"));
        assertFalse(ComplexClassifier.isComplexEvent("APC", "APC p.I1557fs"));

        assertTrue(ComplexClassifier.isComplexEvent("ERBB2", "S310F/Y"));
        assertFalse(ComplexClassifier.isComplexEvent("EGFR", "Exon 19 deletion/insertion"));

        assertFalse(ComplexClassifier.isComplexEvent("BRCA1", "BRCA1 L631QFS"));
    }
}