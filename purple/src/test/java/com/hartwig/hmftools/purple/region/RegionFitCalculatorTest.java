package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.cobalt.CobaltTestUtils;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.purple.purity.RegionFitCalculator;

import org.junit.Test;

public class RegionFitCalculatorTest
{
    private final CobaltChromosomes male = CobaltTestUtils.male();
    private final CobaltChromosomes female = CobaltTestUtils.female();

    @Test
    public void testFitYChromosome()
    {
        final GenomeRegion region = GenomeRegions.create("Y", 1, 100);
        assertTrue(RegionFitCalculator.isAllowedRegion(male, region));
        assertFalse(RegionFitCalculator.isAllowedRegion(female, region));
    }

    @Test
    public void testFitXChromosome()
    {
        final GenomeRegion region = GenomeRegions.create("X", 1, 100);
        assertTrue(RegionFitCalculator.isAllowedRegion(male, region));
        assertTrue(RegionFitCalculator.isAllowedRegion(female, region));
    }
}
