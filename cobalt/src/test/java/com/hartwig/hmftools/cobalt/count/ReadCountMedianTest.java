package com.hartwig.hmftools.cobalt.count;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;

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
    public void testInterpolatedMedian1() {
        var val = Doubles.interpolatedMedian(Arrays.stream((new double[] {  2, 3, 1, 4, 5 })).boxed().collect(Collectors.toList()));
        assertEquals(3.0, val, 1e-10);
    }

    @Test
    public void testModMedian2() {
        addReads(new int[] { 1, 2, 3, 4 });
        assertModMedian(2.5);
    }

    @Test
    public void testInterpolatedMedian2() {
        var val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 3.0, 4.0 })).boxed().collect(Collectors.toList()));
        assertEquals(2.5, val, 1e-10);
    }

    @Test
    public void testModMedian3() {
        addReads(new int[] { 1, 2, 2, 3, 3 });
        assertModMedian(2.25);
    }

    @Test
    public void testInterpolatedMedian3() {
        var val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 2.0, 3.0, 3.0 })).boxed().collect(Collectors.toList()));
        assertEquals(2.25, val, 1e-10);
    }

    @Test
    public void testModMedian4() {
        addReads(new int[] { 1, 2, 2, 2, 3, 3 });
        assertModMedian(2.16666666667);
    }

    @Test
    public void testInterpolatedMedian4() {
        var val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 2.0, 2.0, 3.0, 3.0 })).boxed().collect(Collectors.toList()));
        assertEquals(2.16666666667, val, 1e-10);
    }

    @Test
    public void testModMedian5() {
        addReads(new int[] { 1, 2, 2, 2, 3, 3, 3, 4 });
        assertModMedian(2.5);
    }

    @Test
    public void testInterpolatedMedian5() {
        var val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0 })).boxed().collect(Collectors.toList()));
        assertEquals(2.5, val, 1e-10);
    }
}
