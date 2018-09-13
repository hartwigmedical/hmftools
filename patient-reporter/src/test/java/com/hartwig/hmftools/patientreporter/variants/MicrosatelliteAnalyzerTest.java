package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MicrosatelliteAnalyzerTest {

    @Test
    public void testShortRepeatContextRelevance() {
        assertEquals(false, MicrosatelliteAnalyzer.repeatContextIsRelevant(4, "A"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(5, "A"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(100, "A"));
    }

    @Test
    public void testLongRepeatContextRelevance() {
        assertEquals(false, MicrosatelliteAnalyzer.repeatContextIsRelevant(3, "AT"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(4, "AT"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(100, "AT"));

        assertEquals(false, MicrosatelliteAnalyzer.repeatContextIsRelevant(3, "ATG"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(4, "ATG"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(100, "ATG"));

        assertEquals(false, MicrosatelliteAnalyzer.repeatContextIsRelevant(3, "ATGA"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(4, "ATGA"));
        assertEquals(true, MicrosatelliteAnalyzer.repeatContextIsRelevant(100, "ATGA"));

        assertEquals(false, MicrosatelliteAnalyzer.repeatContextIsRelevant(3, "ATGAA"));
        assertEquals(false, MicrosatelliteAnalyzer.repeatContextIsRelevant(4, "ATGAA"));
    }
}
