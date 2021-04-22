package com.hartwig.hmftools.cobalt.ratio;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class RollingMedianTest {

    private static final double EPSILON = 1e-10;

    private RollingMedian victim;

    @Before
    public void setup() {
        victim = new RollingMedian();
    }

    @Test
    public void testSingleElement() {
        victim.add(5);
        assertMedian(5);
    }

    @Test
    public void testTwoElements() {
        victim.add(5);
        victim.add(6);
        assertMedian(5.5);
    }

    @Test
    public void testThreeElements() {
        victim.add(5);
        victim.add(6);
        victim.add(7);
        assertMedian(6);
    }

    @Test
    public void testRemove() {
        testThreeElements();
        victim.remove(7);
        assertMedian(5.5);
    }

    private void assertMedian(double expected) {
        assertEquals(expected, victim.median(), EPSILON);
    }
}
