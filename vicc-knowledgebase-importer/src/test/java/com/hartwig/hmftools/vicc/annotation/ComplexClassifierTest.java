package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ComplexClassifierTest {

    @Test
    public void canAssessWhetherEventIsComplexEvent() {
        assertTrue(ComplexClassifier.isComplexEvent("Splicing alteration (c.464-2A>G)", "VHL"));

        assertTrue(ComplexClassifier.isComplexEvent("KRAS .", "KRAS"));
        assertTrue(ComplexClassifier.isComplexEvent("APC p.I1557*fs*1", "APC"));
        assertTrue(ComplexClassifier.isComplexEvent("BRCA1 L631Qfs*4", "BRCA1"));
        assertFalse(ComplexClassifier.isComplexEvent("APC p.I1557fs", "APC"));

        assertTrue(ComplexClassifier.isComplexEvent("S310F/Y", "ERBB2"));
        assertFalse(ComplexClassifier.isComplexEvent("Exon 19 deletion/insertion", "EGFR"));

        assertFalse(ComplexClassifier.isComplexEvent("BRCA1 L631QFS", "BRCA1"));
    }
}