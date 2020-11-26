package com.hartwig.hmftools.common.serve.classification.matchers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ComplexMatcherTest {

    @Test
    public void canAssessWhetherEventIsComplexEvent() {
        EventMatcher matcher = new ComplexMatcher();

        assertTrue(matcher.matches("VHL", "Splicing alteration (c.464-2A>G)"));
        assertTrue(matcher.matches("KRAS", "KRAS ."));
        assertTrue(matcher.matches("APC", "APC p.I1557*fs*1"));
        assertTrue(matcher.matches("BRCA1", "BRCA1 L631Qfs*4"));

        assertFalse(matcher.matches("APC", "APC p.I1557fs"));

        assertTrue(matcher.matches("ERBB2", "S310F/Y"));
        assertFalse(matcher.matches("EGFR", "Exon 19 deletion/insertion"));

        assertFalse(matcher.matches("BRCA1", "BRCA1 L631QFS"));
    }
}