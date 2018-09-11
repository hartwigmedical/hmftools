package com.hartwig.hmftools.purple.config;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ConfigSupplierTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testDefaultValues() {
        assertEquals(0.08, FittingConfig.MIN_PURITY_DEFAULT, EPSILON);
        assertEquals(1.0, FittingConfig.MAX_PURITY_DEFAULT, EPSILON);
        assertEquals(0.33, FittingConfig.MIN_NORM_FACTOR_DEFAULT, EPSILON);
        assertEquals(2.0, FittingConfig.MAX_NORM_FACTOR_DEFAULT, EPSILON);
    }
}
