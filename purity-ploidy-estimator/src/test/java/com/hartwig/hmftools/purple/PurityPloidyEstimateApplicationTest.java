package com.hartwig.hmftools.purple;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PurityPloidyEstimateApplicationTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testDefaultValues() {
        assertEquals(0.08, PurityPloidyEstimateApplication.MIN_PURITY_DEFAULT, EPSILON);
        assertEquals(1.0, PurityPloidyEstimateApplication.MAX_PURITY_DEFAULT, EPSILON);
        assertEquals(0.33, PurityPloidyEstimateApplication.MIN_NORM_FACTOR_DEFAULT, EPSILON);
        assertEquals(2.0, PurityPloidyEstimateApplication.MAX_NORM_FACTOR_DEFAULT, EPSILON);
    }
}
