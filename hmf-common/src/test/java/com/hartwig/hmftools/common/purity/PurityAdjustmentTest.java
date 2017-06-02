package com.hartwig.hmftools.common.purity;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedCopyNumber;
import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedVAF;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PurityAdjustmentTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testPurityAdjustedCopynumber() {
        assertEquals(1, purityAdjustedCopyNumber(0.85, 1, 0.575), EPSILON);
    }

    @Test
    public void testPurityAdjustedVAF() {
        System.out.println(purityAdjustedVAF(0.58, 4, 0.34));
        System.out.println(purityAdjustedVAF(0.58, 2, 0.14));
        System.out.println(purityAdjustedVAF(0.58, 2, 0.24));
    }
}
