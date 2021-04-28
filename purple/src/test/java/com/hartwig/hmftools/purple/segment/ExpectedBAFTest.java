package com.hartwig.hmftools.purple.segment;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ExpectedBAFTest
{

    private static final double EPSILON = 1e-8;

    @Test
    public void testExpectedBAF()
    {
        assertEquals(0.56461643, ExpectedBAF.expectedBAF(30, 0.75), EPSILON);
        assertEquals(0.54504357, ExpectedBAF.expectedBAF(60, 0.75), EPSILON);
        assertEquals(0.53681747, ExpectedBAF.expectedBAF(90, 0.75), EPSILON);
        assertEquals(0.53193270, ExpectedBAF.expectedBAF(120, 0.75), EPSILON);
    }
}
