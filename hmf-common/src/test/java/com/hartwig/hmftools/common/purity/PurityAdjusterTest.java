package com.hartwig.hmftools.common.purity;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class PurityAdjusterTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testPurityAdjustedCopynumber() {
        assertPurityAdjustment(0, 0.85, 1.04, 0);
        assertPurityAdjustment(1, 0.85, 1.0, 0.575);
    }

    private void assertPurityAdjustment(final double expectedAdjustedCopyNumber, final double purity, final double normFactor, final double ratio) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, purity, normFactor);
        assertEquals(expectedAdjustedCopyNumber, purityAdjuster.purityAdjustedCopyNumber("1", ratio), EPSILON);
    }

    private PurityAdjuster adjuster(final double purity, final double normFactor) {
        return new PurityAdjuster(Gender.MALE, purity, normFactor);
    }

}
