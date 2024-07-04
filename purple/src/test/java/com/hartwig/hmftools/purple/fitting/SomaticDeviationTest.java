package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.purple.MiscTestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.purple.fittingsnv.SomaticDeviation;

import org.junit.Test;

public class SomaticDeviationTest
{
    @Test
    public void testMaxPloidy()
    {
        PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.FEMALE, 0.12, 0.98);
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
        PurityAdjuster purityAdjuster = buildPurityAdjuster(Gender.FEMALE, 0.12, 0.98);
        int maxReads = SomaticDeviation.INSTANCE.maxConceivableReads(purityAdjuster, 2, depth(18, 55), -0.1, -0.1);
        assertEquals(0, maxReads);
    }

    private static AllelicDepth depth(int alleleReadCount, int totalReadCount)
    {
        return new AllelicDepth(totalReadCount, alleleReadCount);
    }
}
