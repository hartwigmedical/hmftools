package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ComplexClassifierTest {

    @Test
    public void canAssessWhetherEventIsComplexEvent() {
        assertTrue(ComplexClassifier.isComplexEvent("B2M .", "B2M"));
        assertTrue(ComplexClassifier.isComplexEvent("APC p.I1557*fs*1", "APC"));
        assertTrue(ComplexClassifier.isComplexEvent("BRCA1 L631Qfs*4", "BRCA1"));
    }

}