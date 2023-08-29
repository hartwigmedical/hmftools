package com.hartwig.hmftools.purple.purity;

import static com.hartwig.hmftools.purple.TestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.config.FittingConfig;

import org.junit.Test;

public class PurityAdjusterTest
{
    private static final double EPSILON = 1e-10;

    @Test
    public void testPurityAdjustedCopyNumber()
    {
        assertPurityAdjustment(0, 0.85, 1.04, 0);
        assertPurityAdjustment(1, 0.85, 1.0, 0.575);
    }

    private static void assertPurityAdjustment(final double expectedAdjustedCopyNumber, final double purity, final double normFactor,
            final double ratio)
    {
        final PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.MALE, purity, normFactor);
        assertEquals(expectedAdjustedCopyNumber, purityAdjuster.purityAdjustedCopyNumber("1", ratio), EPSILON);
    }

    @Test
    public void testPurityAdjustedFrequency()
    {
        final PurityAdjuster victim = buildPurityAdjuster(Gender.FEMALE, 0.9, 1d);

        victim.purityAdjustedFrequency(2, 1, 3.0, 0.33);

        assertFrequencyMatchesPloidy(victim, 2, 1, 3, 2);
        assertFrequencyMatchesPloidy(victim, 2, 0, 3, 2);
    }

    @Test
    public void testPureTumor()
    {
        final PurityAdjuster victim = buildPurityAdjuster(Gender.FEMALE, 1, 1d);
        assertFrequencyMatchesPloidy(victim, 2, 1, 3, 2);
        assertFrequencyMatchesPloidy(victim, 2, 0, 3, 2);
    }

    @Test
    public void testExpectedFrequencyWithNegativeTumorCopyNumber()
    {
        final PurityAdjuster victim = buildPurityAdjuster(Gender.FEMALE, 1, 1d);
        assertEquals(0, victim.expectedFrequency(2, 1, -0.2, -0.1), EPSILON);
    }

    @Test
    public void testNormFactorPloidyAverageRatioConsistency()
    {
        double averageRatio = 0.990374738586366;
        double normFactor = 0.8782;
        double purity = 0.08;

        double impliedPloidy = PurityAdjuster.impliedPloidy(averageRatio, purity, normFactor);
        double impliedNormFactor = PurityAdjuster.impliedNormFactor(averageRatio, purity, impliedPloidy);
        double impliedAverageRatio = PurityAdjuster.impliedAverageRatio(purity, normFactor, impliedPloidy);

        assertEquals(averageRatio, impliedAverageRatio, EPSILON);
        assertEquals(normFactor, impliedNormFactor, EPSILON);
    }

    @Test
    public void testBafDeviationCalcs()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        FittingConfig.addConfig(configBuilder);
        FittingConfig fittingConfig = new FittingConfig(configBuilder);

        String chromosome = "1";
        MedianRatio medianRatio = new MedianRatio(chromosome, 0.5, 1);
        CobaltChromosomes cobaltChromosomes = new CobaltChromosomes(Lists.newArrayList(medianRatio), true);

        RegionFitCalculator regionFitCalculator = new RegionFitCalculator(cobaltChromosomes, fittingConfig, 100);
        double purity = 0.08;
        double normFactor = 0.986;

        PurityAdjuster purityAdjuster = new PurityAdjuster(purity, normFactor, cobaltChromosomes);

        double baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 2.018, 0.53);
        assertEquals(0.504, baf, 0.001);

        baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 1.5, 0.53);
        assertEquals(0.667, baf, 0.001);

        baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 2.02, 0.51);
        assertEquals(0.504, baf, 0.001);

        baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 1.9, 0.53);
        assertEquals(0.526, baf, 0.001);
    }

    private static void assertFrequencyMatchesPloidy(
            final PurityAdjuster victim, final int normalCopyNumber,
            final int normalPloidy, final double tumorCopyNumber, final double tumorPloidy)
    {
        double expectedFrequency = victim.expectedFrequency(normalCopyNumber, normalPloidy, tumorCopyNumber, tumorPloidy);
        double actualTumorPloidy = victim.purityAdjustedPloidy(normalCopyNumber, normalPloidy, tumorCopyNumber, expectedFrequency);

        assertEquals(tumorPloidy, actualTumorPloidy, EPSILON);
    }
}