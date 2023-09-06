package com.hartwig.hmftools.common.genome.chromosome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class MitochondrialChromosomeTest
{
    @Test
    public void testFromString()
    {
        assertEquals(MitochondrialChromosome.MT, MitochondrialChromosome.fromString("MT"));
        assertEquals(MitochondrialChromosome.MT, MitochondrialChromosome.fromString("chrM"));
    }

    @Test
    public void testContains()
    {
        assertTrue(MitochondrialChromosome.contains("MT"));
        assertTrue(MitochondrialChromosome.contains("chrM"));

        assertFalse(MitochondrialChromosome.contains("1"));
        assertFalse(MitochondrialChromosome.contains("chr1"));
        assertFalse(MitochondrialChromosome.contains("X"));
        assertFalse(MitochondrialChromosome.contains("chrX"));
    }
}
