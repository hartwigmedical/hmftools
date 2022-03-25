package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.cobalt.CobaltTestUtils;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.junit.Test;

public class FittedRegionFactoryTest
{
    private final CobaltChromosomes male = CobaltTestUtils.male();
    private final CobaltChromosomes female = CobaltTestUtils.female();

    @Test
    public void testFitYChromosome()
    {
        final GenomeRegion region = GenomeRegions.create("Y", 1, 100);
        assertTrue(FittedRegionFactory.isAllowedRegion(male, region));
        assertFalse(FittedRegionFactory.isAllowedRegion(female, region));
    }

    @Test
    public void testFitXChromosome()
    {
        final GenomeRegion region = GenomeRegions.create("X", 1, 100);
        assertTrue(FittedRegionFactory.isAllowedRegion(male, region));
        assertTrue(FittedRegionFactory.isAllowedRegion(female, region));
    }
}
