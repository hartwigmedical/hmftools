package com.hartwig.hmftools.common.purple.purity;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;

import org.junit.Test;

public class SomaticDeviationTest {

    @Test
    public void testMaxPloidy() {
        PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.FEMALE, 0.12, 0.98);
        int maxReads = SomaticDeviation.INSTANCE.maxConceivableReads(purityAdjuster, 2, depth(18, 55),  3.0965, 2.2021);
        assertEquals(14, maxReads);

        double maxPloidy = SomaticDeviation.INSTANCE.maxConceivablePloidy(purityAdjuster, 2, depth(18, 55),  3.0965, 2.2021);
        assertEquals(4.52, maxPloidy, 0.01);

        double deviationFromMax = SomaticDeviation.INSTANCE.deviationFromMax(purityAdjuster, 2, depth(18, 55),  3.0965, 2.2021);
        assertEquals(1.29, deviationFromMax, 0.01);
    }

    private AllelicDepth depth(int alleleReadCount, int totalReadCount) {
        return ImmutableAllelicDepthImpl.builder().alleleReadCount(alleleReadCount).totalReadCount(totalReadCount).build();
    }


}
