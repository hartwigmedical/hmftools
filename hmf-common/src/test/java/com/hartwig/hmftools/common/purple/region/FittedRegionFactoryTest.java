package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class FittedRegionFactoryTest {

    @Test
    public void testFitYChromosome() {
        final GenomeRegion region = GenomeRegions.create("Y", 1, 100);
        assertTrue(FittedRegionFactoryV2.isFittableRegion(Gender.MALE, region));
        assertTrue(FittedRegionFactoryV2.isFittableRegion(Gender.MALE_KLINEFELTER, region));
        assertFalse(FittedRegionFactoryV2.isFittableRegion(Gender.FEMALE, region));
    }

    @Test
    public void testFitChrYChromosome() {
        final GenomeRegion region = GenomeRegions.create("chrY", 1, 100);
        assertTrue(FittedRegionFactoryV2.isFittableRegion(Gender.MALE, region));
        assertTrue(FittedRegionFactoryV2.isFittableRegion(Gender.MALE_KLINEFELTER, region));
        assertFalse(FittedRegionFactoryV2.isFittableRegion(Gender.FEMALE, region));
    }

    @Test
    public void testFitXChromosome() {
        final GenomeRegion region = GenomeRegions.create("X", 1, 100);
        assertTrue(FittedRegionFactoryV2.isFittableRegion(Gender.MALE, region));
        assertTrue(FittedRegionFactoryV2.isFittableRegion(Gender.MALE_KLINEFELTER, region));
        assertTrue(FittedRegionFactoryV2.isFittableRegion(Gender.FEMALE, region));
    }

}
