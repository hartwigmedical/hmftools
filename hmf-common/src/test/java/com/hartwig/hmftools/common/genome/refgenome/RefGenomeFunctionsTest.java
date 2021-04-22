package com.hartwig.hmftools.common.genome.refgenome;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RefGenomeFunctionsTest {

    @Test
    public void canStripChromosomes() {
        assertEquals("10", RefGenomeFunctions.stripChromosome("chr10"));
        assertEquals("10", RefGenomeFunctions.stripChromosome("10"));
    }

    @Test
    public void canEnforceChromosomes() {
        assertEquals("chr10", RefGenomeFunctions.enforceChromosome("chr10"));
        assertEquals("chr10", RefGenomeFunctions.enforceChromosome("10"));
    }
}