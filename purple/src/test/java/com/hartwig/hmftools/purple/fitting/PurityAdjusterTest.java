package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.FittingTestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.Gender;

import org.junit.Test;

public class PurityAdjusterTest
{
    private static final double EPSILON = 1e-10;

    @Test
    public void germlineCopyNumber()
    {
        final PurityAdjuster adjuster = buildPurityAdjuster(Gender.FEMALE, 1, 1d);
        assertEquals(2.0, adjuster.germlineCopyNumber("1"), EPSILON);
        assertEquals(2.0, adjuster.germlineCopyNumber("X"), EPSILON);
        assertEquals(0.0, adjuster.germlineCopyNumber("Y"), EPSILON);
    }

    @Test
    public void testPurityAdjustedCopyNumber()
    {
        PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.MALE, 1.0, 1.0);
        assertEquals(2.0, purityAdjuster.purityAdjustedCopyNumber("1", 1.0), EPSILON);
        assertEquals(2.4, purityAdjuster.purityAdjustedCopyNumber("1", 1.2), EPSILON);
        assertEquals(4.0, purityAdjuster.purityAdjustedCopyNumber("1", 2.0), EPSILON);

        purityAdjuster = buildPurityAdjuster(Gender.MALE, 0.80, 1.0);
        assertEquals(2.0, purityAdjuster.purityAdjustedCopyNumber("1", 1.0), EPSILON);
        assertEquals(3.25, purityAdjuster.purityAdjustedCopyNumber("1", 1.5), EPSILON);
        assertEquals(4.5, purityAdjuster.purityAdjustedCopyNumber("1", 2.0), EPSILON);

        purityAdjuster = buildPurityAdjuster(Gender.MALE, 0.85, 1.04);
        assertEquals(0.0, purityAdjuster.purityAdjustedCopyNumber("1", 0.0), EPSILON);

        purityAdjuster = buildPurityAdjuster(Gender.MALE, 0.85, 1.00);
        assertEquals(1.0, purityAdjuster.purityAdjustedCopyNumber("1", 0.575), EPSILON);
    }

    @Test
    public void testPurityAdjustedCopyNumberGivenRatio()
    {
        // 2 * normalRatio + 2 * (tumorRatio - normalRatio * mNormFactor) / mPurity / mNormFactor;
        PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.MALE, 1.0, 1.0);
        assertEquals(2.0, purityAdjuster.purityAdjustedCopyNumber(1.0, 1.0), EPSILON);
        assertEquals(4.0, purityAdjuster.purityAdjustedCopyNumber(2.0, 1.0), EPSILON);
        assertEquals(2.0, purityAdjuster.purityAdjustedCopyNumber(1.0, 2.0), EPSILON);

        purityAdjuster = buildPurityAdjuster(Gender.MALE, 1.0, 0.5);
        assertEquals(4.0, purityAdjuster.purityAdjustedCopyNumber(1.0, 1.0), EPSILON);
        assertEquals(8.0, purityAdjuster.purityAdjustedCopyNumber(2.0, 1.0), EPSILON);
        assertEquals(4.0, purityAdjuster.purityAdjustedCopyNumber(1.0, 2.0), EPSILON);

        purityAdjuster = buildPurityAdjuster(Gender.MALE, 0.5, 1.0);
        assertEquals(2.0, purityAdjuster.purityAdjustedCopyNumber(1.0, 1.0), EPSILON);
        assertEquals(6.0, purityAdjuster.purityAdjustedCopyNumber(2.0, 1.0), EPSILON);
        assertEquals(0.0, purityAdjuster.purityAdjustedCopyNumber(1.0, 2.0), EPSILON);
        assertEquals(-4.0, purityAdjuster.purityAdjustedCopyNumber(1.0, 4.0), EPSILON);
    }

    @Test
    public void testPurityAdjustedFrequency()
    {
        PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.FEMALE, 0.9, 1d);

        purityAdjuster.purityAdjustedFrequency(2, 1, 3.0, 0.33);

        assertFrequencyMatchesPloidy(purityAdjuster, 2, 1, 3, 2);
        assertFrequencyMatchesPloidy(purityAdjuster, 2, 0, 3, 2);
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
        PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.FEMALE, 1, 1d);
        assertEquals(0, purityAdjuster.expectedFrequency(2, 1, -0.2, -0.1), EPSILON);
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