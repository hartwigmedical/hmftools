package com.hartwig.hmftools.purple.purity;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SomaticDeviationTest
{
    @Test
    public void testMaxPloidy()
    {
        PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.FEMALE, 0.12, 0.98);
        int maxReads = SomaticDeviation.INSTANCE.maxConceivableReads(purityAdjuster, 2, depth(18, 55), 3.0965, 2.2021);
        assertEquals(16, maxReads);

        double maxPloidy = SomaticDeviation.INSTANCE.maxConceivablePloidy(purityAdjuster, 2, depth(18, 55), 3.0965, 2.2021);
        assertEquals(5.16, maxPloidy, 0.01);

        double deviationFromMax = SomaticDeviation.INSTANCE.deviationFromMax(purityAdjuster, 2, depth(18, 55), 3.0965, 2.2021);
        assertEquals(0.65, deviationFromMax, 0.01);
    }

    @Test
    public void testMaxConceivableReadsIsZeroWithNegativeTumorCopyNumber()
    {
        PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.FEMALE, 0.12, 0.98);
        int maxReads = SomaticDeviation.INSTANCE.maxConceivableReads(purityAdjuster, 2, depth(18, 55), -0.1, -0.1);
        assertEquals(0, maxReads);
    }

    @NotNull
    private static AllelicDepth depth(int alleleReadCount, int totalReadCount)
    {
        return ImmutableAllelicDepthImpl.builder().alleleReadCount(alleleReadCount).totalReadCount(totalReadCount).build();
    }
}
