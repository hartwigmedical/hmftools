package com.hartwig.hmftools.cobalt.count;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

public class ReadCountMedianTest {

    private ReadCountMedian victim;

    @Before
    public void setup() {
        victim = new ReadCountMedian();
    }

    @Test
    public void testNoReads() {
        assertMedian(0);
    }

    @Test
    public void testOneRead() {
        victim.addRead(100);
        assertMedian(100);
    }

    @Test
    public void testTwoReads() {
        victim.addRead(100);
        victim.addRead(150);
        assertMedian(125);
    }

    @Test
    public void testThreeReads() {
        victim.addRead(100);
        victim.addRead(150);
        victim.addRead(130);
        assertMedian(130);
    }

    @Test
    public void testFourReads() {
        victim.addRead(100);
        victim.addRead(100);
        victim.addRead(150);
        victim.addRead(130);
        assertMedian(115);
    }

    @Test
    public void testFiveReads() {
        victim.addRead(100);
        victim.addRead(100);
        victim.addRead(150);
        victim.addRead(130);
        victim.addRead(130);
        assertMedian(130);
    }

    private void addReads(int[] reads) {
        for (int i : reads)
            victim.addRead(i);
    }

    private void assertMedian(int expectedMedian) {
        assertEquals(expectedMedian, victim.median(), 1e-10);
    }

    private void assertModMedian(double expectedModMedian) {
        assertEquals(expectedModMedian, victim.interpolatedMedian(), 1e-10);
    }

    @Test
    public void testModMedian1() {
        addReads(new int[] { 2, 3, 1, 4, 5 });
        assertModMedian(3);
    }

    @Test
    public void testModMedian2() {
        addReads(new int[] { 1, 2, 3, 4 });
        assertModMedian(2.5);
    }

    @Test
    public void testModMedian3() {
        addReads(new int[] { 1, 2, 2, 3, 3 });
        assertModMedian(2.25);
    }

    @Test
    public void testModMedian4() {
        addReads(new int[] { 1, 2, 2, 2, 3, 3 });
        assertModMedian(2.16666666667);
    }
}
