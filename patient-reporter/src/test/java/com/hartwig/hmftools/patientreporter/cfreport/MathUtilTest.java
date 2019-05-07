package com.hartwig.hmftools.patientreporter.cfreport;

import static org.junit.Assert.*;

import org.junit.Test;

public class MathUtilTest {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void mapPercentageWorksAsExpected() {
        assertEquals(50D, MathUtil.mapPercentage(0.25, 0, 0.5), EPSILON);
        assertEquals(25D, MathUtil.mapPercentage(1.5, 1, 3), EPSILON);
        assertEquals(-100D, MathUtil.mapPercentage(-1, 1, 3), EPSILON);
    }

    @Test
    public void mapPercentageClampedWorksAsExpected() {
        assertEquals(50D, MathUtil.mapPercentageClamped(0.25, 0, 0.5), EPSILON);
        assertEquals(25D, MathUtil.mapPercentageClamped(1.5, 1, 3), EPSILON);
        assertEquals(0D, MathUtil.mapPercentageClamped(-1, 1, 3), EPSILON);
    }

}