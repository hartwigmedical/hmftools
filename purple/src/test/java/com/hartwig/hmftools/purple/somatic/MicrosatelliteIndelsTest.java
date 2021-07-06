package com.hartwig.hmftools.purple.somatic;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class MicrosatelliteIndelsTest {

    @Test
    public void shortRepeatContextRelevanceIsDeterminedCorrectly() {
        assertFalse(MicrosatelliteIndels.repeatContextIsRelevant(4, "A"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(5, "A"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(100, "A"));
    }

    @Test
    public void longRepeatContextRelevanceIsDeterminedCorrectly() {
        assertFalse(MicrosatelliteIndels.repeatContextIsRelevant(3, "AT"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(4, "AT"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(100, "AT"));

        assertFalse(MicrosatelliteIndels.repeatContextIsRelevant(3, "ATG"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(4, "ATG"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(100, "ATG"));

        assertFalse(MicrosatelliteIndels.repeatContextIsRelevant(3, "ATGA"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(4, "ATGA"));
        assertTrue(MicrosatelliteIndels.repeatContextIsRelevant(100, "ATGA"));

        assertFalse(MicrosatelliteIndels.repeatContextIsRelevant(3, "ATGAA"));
        assertFalse(MicrosatelliteIndels.repeatContextIsRelevant(4, "ATGAA"));
    }
}
