package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;

import org.junit.Test;

public class ComplexClassifierTest {

    @Test
    public void canAssessWhetherEventIsComplexEvent() {
        EventClassifier classifier = new ComplexClassifier();

        assertTrue(classifier.matches("VHL", "Splicing alteration (c.464-2A>G)"));
        assertTrue(classifier.matches("KRAS", "KRAS ."));
        assertTrue(classifier.matches("APC", "APC p.I1557*fs*1"));
        assertTrue(classifier.matches("BRCA1", "BRCA1 L631Qfs*4"));

        assertFalse(classifier.matches("APC", "APC p.I1557fs"));

        assertTrue(classifier.matches("ERBB2", "S310F/Y"));
        assertFalse(classifier.matches("EGFR", "Exon 19 deletion/insertion"));

        assertFalse(classifier.matches("BRCA1", "BRCA1 L631QFS"));
    }
}