package com.hartwig.hmftools.common.msi;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MicrosatelliteIndelsTest {

    @Test
    public void shortRepeatContextRelevanceIsDeterminedCorrectly() {
        assertEquals(false, MicrosatelliteIndels.repeatContextIsRelevant(4, "A"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(5, "A"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(100, "A"));
    }

    @Test
    public void longRepeatContextRelevanceIsDeterminedCorrectly() {
        assertEquals(false, MicrosatelliteIndels.repeatContextIsRelevant(3, "AT"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(4, "AT"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(100, "AT"));

        assertEquals(false, MicrosatelliteIndels.repeatContextIsRelevant(3, "ATG"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(4, "ATG"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(100, "ATG"));

        assertEquals(false, MicrosatelliteIndels.repeatContextIsRelevant(3, "ATGA"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(4, "ATGA"));
        assertEquals(true, MicrosatelliteIndels.repeatContextIsRelevant(100, "ATGA"));

        assertEquals(false, MicrosatelliteIndels.repeatContextIsRelevant(3, "ATGAA"));
        assertEquals(false, MicrosatelliteIndels.repeatContextIsRelevant(4, "ATGAA"));
    }
}
