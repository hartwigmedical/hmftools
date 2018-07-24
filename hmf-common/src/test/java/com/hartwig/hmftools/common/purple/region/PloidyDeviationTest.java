package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PloidyDeviationTest {

    private static final double EPSILON = 1e-3;

    @Test
    public void testLowPurity() {
        final PloidyDeviation victim = new PloidyDeviation(0.4, 0.63, 0.03);
        assertEquals(0.326, victim.deviation(0.1), EPSILON);
        assertEquals(0.599, victim.deviation(0.2), EPSILON);
        assertEquals(0.792, victim.deviation(0.3), EPSILON);
        assertEquals(0.907, victim.deviation(0.4), EPSILON);
    }

    @Test
    public void testHighPurity() {
        final PloidyDeviation victim = new PloidyDeviation(1, 0.63, 0.03);
        assertEquals(0.706, victim.deviation(0.1), EPSILON);
        assertEquals(0.964, victim.deviation(0.2), EPSILON);
        assertEquals(0.998, victim.deviation(0.3), EPSILON);
        assertEquals(1.000, victim.deviation(0.4), EPSILON);
    }

}
