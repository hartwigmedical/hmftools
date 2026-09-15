package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._7;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class CnBackboneResolutionTest
{
    @Test
    public void testNoOverrides()
    {
        CnBackboneResolution resolution = CnBackboneResolution.parse(1000, null);
        assertEquals(1_000_000, resolution.forChromosome(_7));
        assertEquals(1_000_000, resolution.forChromosome(_Y));
    }

    // Overridden chromosomes use their own spacing; others fall through to the default.
    @Test
    public void testOverrides()
    {
        CnBackboneResolution resolution = CnBackboneResolution.parse(1000, "7:500,Y:2000");
        assertEquals(500_000, resolution.forChromosome(_7));
        assertEquals(2_000_000, resolution.forChromosome(_Y));
        assertEquals(1_000_000, resolution.forChromosome(HumanChromosome._1));
    }

    @Test
    public void testChrPrefixAccepted()
    {
        CnBackboneResolution resolution = CnBackboneResolution.parse(1000, "chr7:500");
        assertEquals(500_000, resolution.forChromosome(_7));
    }

    @Test
    public void testDefaultBelowMinimumRejected()
    {
        assertThrows(UserInputError.class, () -> CnBackboneResolution.parse(0, null));
    }

    @Test
    public void testOverrideBelowMinimumRejected()
    {
        assertThrows(UserInputError.class, () -> CnBackboneResolution.parse(1000, "7:0"));
    }

    @Test
    public void testUnknownChromosomeRejected()
    {
        assertThrows(UserInputError.class, () -> CnBackboneResolution.parse(1000, "25:500"));
    }

    @Test
    public void testDuplicateChromosomeRejected()
    {
        assertThrows(UserInputError.class, () -> CnBackboneResolution.parse(1000, "7:500,7:600"));
    }

    @Test
    public void testMalformedPairRejected()
    {
        assertThrows(UserInputError.class, () -> CnBackboneResolution.parse(1000, "7-500"));
    }

    @Test
    public void testNonNumericResolutionRejected()
    {
        assertThrows(UserInputError.class, () -> CnBackboneResolution.parse(1000, "7:abc"));
    }
}
