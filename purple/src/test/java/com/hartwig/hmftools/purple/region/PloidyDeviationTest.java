package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PloidyDeviationTest
{

    private static final double EPSILON = 1e-3;

    @Test
    public void testLowPurity()
    {
        final PloidyDeviation victim = new PloidyDeviation(0.03, 0, 1, 1, 0);
        assertEquals(0.326, victim.minorAlleleDeviation(0.4, 0.63, 0.1), EPSILON);
        assertEquals(0.599, victim.minorAlleleDeviation(0.4, 0.63, 0.2), EPSILON);
        assertEquals(0.792, victim.minorAlleleDeviation(0.4, 0.63, 0.3), EPSILON);
        assertEquals(0.907, victim.minorAlleleDeviation(0.4, 0.63, 0.4), EPSILON);
        assertEquals(0.964, victim.minorAlleleDeviation(0.4, 0.63, 0.5), EPSILON);
        assertEquals(0.907, victim.minorAlleleDeviation(0.4, 0.63, 0.6), EPSILON);
        assertEquals(0.792, victim.minorAlleleDeviation(0.4, 0.63, 0.7), EPSILON);
        assertEquals(0.599, victim.minorAlleleDeviation(0.4, 0.63, 0.8), EPSILON);
        assertEquals(0.326, victim.minorAlleleDeviation(0.4, 0.63, 0.9), EPSILON);
        assertEquals(0.000, victim.majorAlleleDeviation(0.4, 0.63, 1.0), EPSILON);
        assertEquals(0.326, victim.minorAlleleDeviation(0.4, 0.63, 1.1), EPSILON);
    }

    @Test
    public void testMinDeviation()
    {
        double minDeviation = 0.4;
        final PloidyDeviation victim = new PloidyDeviation(0.03, 0, 1, 1, minDeviation);
        assertEquals(Math.max(minDeviation, 0.326), victim.minorAlleleDeviation(0.4, 0.63, 0.1), EPSILON);
        assertEquals(Math.max(minDeviation, 0.599), victim.minorAlleleDeviation(0.4, 0.63, 0.2), EPSILON);
        assertEquals(Math.max(minDeviation, 0.792), victim.minorAlleleDeviation(0.4, 0.63, 0.3), EPSILON);
        assertEquals(Math.max(minDeviation, 0.907), victim.minorAlleleDeviation(0.4, 0.63, 0.4), EPSILON);
        assertEquals(Math.max(minDeviation, 0.964), victim.minorAlleleDeviation(0.4, 0.63, 0.5), EPSILON);
        assertEquals(Math.max(minDeviation, 0.907), victim.minorAlleleDeviation(0.4, 0.63, 0.6), EPSILON);
        assertEquals(Math.max(minDeviation, 0.792), victim.minorAlleleDeviation(0.4, 0.63, 0.7), EPSILON);
        assertEquals(Math.max(minDeviation, 0.599), victim.minorAlleleDeviation(0.4, 0.63, 0.8), EPSILON);
        assertEquals(Math.max(minDeviation, 0.326), victim.minorAlleleDeviation(0.4, 0.63, 0.9), EPSILON);
        assertEquals(Math.max(minDeviation, 0.000), victim.majorAlleleDeviation(0.4, 0.63, 1.0), EPSILON);
        assertEquals(Math.max(minDeviation, 0.326), victim.minorAlleleDeviation(0.4, 0.63, 1.1), EPSILON);
    }

    @Test
    public void testHighPurity()
    {
        final PloidyDeviation victim = new PloidyDeviation(0.03, 0, 1, 1, 0);
        assertEquals(0.706, victim.minorAlleleDeviation(1, 0.63, 0.1), EPSILON);
        assertEquals(0.964, victim.minorAlleleDeviation(1, 0.63, 0.2), EPSILON);
        assertEquals(0.998, victim.minorAlleleDeviation(1, 0.63, 0.3), EPSILON);
        assertEquals(1.000, victim.minorAlleleDeviation(1, 0.63, 0.4), EPSILON);
    }
}
