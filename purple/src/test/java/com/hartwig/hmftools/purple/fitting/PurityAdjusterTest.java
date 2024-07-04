package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.MiscTestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.Gender;

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

    private static void assertFrequencyMatchesPloidy(
            final PurityAdjuster victim, final int normalCopyNumber,
            final int normalPloidy, final double tumorCopyNumber, final double tumorPloidy)
    {
        double expectedFrequency = victim.expectedFrequency(normalCopyNumber, normalPloidy, tumorCopyNumber, tumorPloidy);
        double actualTumorPloidy = victim.purityAdjustedPloidy(normalCopyNumber, normalPloidy, tumorCopyNumber, expectedFrequency);

        assertEquals(tumorPloidy, actualTumorPloidy, EPSILON);
    }
}