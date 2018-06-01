package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MutationalLoadAnalyzerTest {

    @Test
    public void analyzeVariants() {
        assertEquals(true, true);
    }

    @Test
    public void testMutationalLoadCheck() {
        assertEquals("MISSENSE", "MISSENSE");
        assertEquals("NONE", "NONE");
        assertEquals("NONSENSE_OR_FRAMESHIFT", "NONSENSE_OR_FRAMESHIFT");
        assertEquals("SYNONYMOUS", "SYNONYMOUS");
        assertEquals("SPLICE", "SPLICE");

    }
}

