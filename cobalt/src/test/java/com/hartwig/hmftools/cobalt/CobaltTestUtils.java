package com.hartwig.hmftools.cobalt;

import static org.junit.Assert.assertEquals;

public class CobaltTestUtils
{
    public static final double EPSILON = 1e-7;

    public static void assertDoubleEquals(double expected, double actual)
    {
        assertEquals(expected, actual, EPSILON);
    }
}
