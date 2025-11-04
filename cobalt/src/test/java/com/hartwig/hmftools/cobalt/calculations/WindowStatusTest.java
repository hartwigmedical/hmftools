package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class WindowStatusTest extends CalculationsTestBase
{

    @Test
    public void maskedOutByMappability()
    {
        assertFalse(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.99), false, false).maskedOut());
        assertFalse(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.87), false, false).maskedOut());
        assertFalse(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.86), false, false).maskedOut());
        assertFalse(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.85), false, false).maskedOut());
        assertFalse(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.85), false, false).maskedOut());
        assertTrue(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.84), false, false).maskedOut());
        assertTrue(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.83), false, false).maskedOut());
        assertTrue(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.01), false, false).maskedOut());
    }

    @Test
    public void maskedOutByExlusion()
    {
        assertTrue(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.99), true, false).maskedOut());
        assertFalse(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.87), false, false).maskedOut());
        assertTrue(new WindowStatus(gcProfile(_1, 1001, 0.47, 0.86), true, false).maskedOut());
    }
}
