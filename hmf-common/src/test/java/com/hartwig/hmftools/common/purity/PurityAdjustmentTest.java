package com.hartwig.hmftools.common.purity;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedBAF;
import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedCopynumber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;

import org.junit.Test;

public class PurityAdjustmentTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testPurityAdjustedCopynumber() {
        assertEquals(1, purityAdjustedCopynumber(0.85, 1, 0.575), EPSILON);
    }

    @Test
    public void testPurityAdjustedBaf() {
        testPurityAdjustedBaf(0.8, 3, 1);
        testPurityAdjustedBaf(0.8, 3, 2);
        testPurityAdjustedBaf(0.9, 4, 1);
    }

    private void testPurityAdjustedBaf(double purity, int ploidy, int alleleCount) {
        double expectedPurityAdjustedBAF = 1d * alleleCount / ploidy;
        double observedBAF = FittedRegionFactory.modelBAF(purity, ploidy, alleleCount);
        assertEquals(expectedPurityAdjustedBAF, purityAdjustedBAF(purity, ploidy, observedBAF), EPSILON);
    }
}
