package com.hartwig.hmftools.common.purity;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedCopyNumber;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PurityAdjustmentTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testPurityAdjustedCopynumber() {
        assertEquals(0, purityAdjustedCopyNumber(0.85, 1.04, 0, 2), EPSILON);
        assertEquals(1, purityAdjustedCopyNumber(0.85, 1, 0.575, 2), EPSILON);
    }
}
