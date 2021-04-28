package com.hartwig.hmftools.purple.config;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PurpleConfigTest
{

    private static final double EPSILON = 1e-10;

    @Test
    public void testDefaultValues()
    {
        assertEquals(0.08, FittingConfig.MIN_PURITY_DEFAULT, EPSILON);
        assertEquals(1.0, FittingConfig.MAX_PURITY_DEFAULT, EPSILON);
        assertEquals(1, FittingConfig.MIN_PLOIDY_DEFAULT, EPSILON);
        assertEquals(8, FittingConfig.MAX_PLOIDY_DEFAULT, EPSILON);
    }
}
