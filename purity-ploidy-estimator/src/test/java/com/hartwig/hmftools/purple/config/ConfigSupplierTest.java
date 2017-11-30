package com.hartwig.hmftools.purple.config;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ConfigSupplierTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testDefaultValues() {
        assertEquals(0.08, ConfigSupplier.MIN_PURITY_DEFAULT, EPSILON);
        assertEquals(1.0, ConfigSupplier.MAX_PURITY_DEFAULT, EPSILON);
        assertEquals(0.33, ConfigSupplier.MIN_NORM_FACTOR_DEFAULT, EPSILON);
        assertEquals(2.0, ConfigSupplier.MAX_NORM_FACTOR_DEFAULT, EPSILON);
    }
}
