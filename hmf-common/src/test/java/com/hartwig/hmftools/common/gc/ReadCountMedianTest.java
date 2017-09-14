package com.hartwig.hmftools.common.gc;

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

    private void assertMedian(int expectedMedian) {
        assertEquals(expectedMedian, victim.median());
    }

}
