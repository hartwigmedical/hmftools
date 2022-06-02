package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ChromosomesTest {

    @Test
    public void canPrefixChromosomes() {
        assertEquals("05", Chromosomes.zeroPrefixed("5"));
        assertEquals("05p13.2", Chromosomes.zeroPrefixed("5p13.2"));

        assertEquals("15", Chromosomes.zeroPrefixed("15"));
        assertEquals("15q21.1", Chromosomes.zeroPrefixed("15q21.1"));

        assertEquals("X", Chromosomes.zeroPrefixed("X"));
    }
}