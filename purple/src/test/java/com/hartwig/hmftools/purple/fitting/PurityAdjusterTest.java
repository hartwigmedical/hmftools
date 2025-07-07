package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.FittingTestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;

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
        PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.MALE, purity, normFactor);
        assertEquals(expectedAdjustedCopyNumber, purityAdjuster.purityAdjustedCopyNumber("1", ratio), EPSILON);
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

    @Test
    public void testMostDiploidPurity()
    {
        FittedPurity fp1 = createFittedPurity(0.3, 0.3, 2.3);
        FittedPurity fp2 = createFittedPurity(0.3, 0.2, 1.9);
        FittedPurity fp3 = createFittedPurity(0.4, 0.4, 1.8);
        FittedPurity fp4 = createFittedPurity(0.4, 0.3, 2.05);

        List<FittedPurity> all = Lists.newArrayList(fp1, fp2, fp3, fp4);
        Collections.shuffle(all);

        List<FittedPurity> result = BestFit.mostDiploidPerPurity(all);

        assertEquals(2, result.size());
        assertEquals(fp2, result.get(0));
        assertEquals(fp4, result.get(1));
    }

    private FittedPurity createFittedPurity(double purity, double score, double ploidy)
    {
        return ImmutableFittedPurity.builder()
                .purity(purity)
                .normFactor(0.95)
                .score(score)
                .diploidProportion(1.0)
                .ploidy(ploidy)
                .somaticPenalty(0).build();
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