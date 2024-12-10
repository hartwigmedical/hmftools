package com.hartwig.hmftools.cobalt.count;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;

import org.junit.Test;

public class ReadCountMedianTest
{
    @Test
    public void testInterpolatedMedian1()
    {
        double val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 2, 3, 1, 4, 5 })).boxed().collect(Collectors.toList()));
        assertEquals(3.0, val, 1e-10);
    }

    @Test
    public void testInterpolatedMedian2()
    {
        double val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 3.0, 4.0 })).boxed().collect(Collectors.toList()));
        assertEquals(2.5, val, 1e-10);
    }

    @Test
    public void testInterpolatedMedian3()
    {
        double val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 2.0, 3.0, 3.0 })).boxed().collect(Collectors.toList()));
        assertEquals(2.25, val, 1e-10);
    }

    @Test
    public void testInterpolatedMedian4()
    {
        double val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 2.0, 2.0, 3.0, 3.0 }))
                .boxed()
                .collect(Collectors.toList()));
        assertEquals(2.16666666667, val, 1e-10);
    }

    @Test
    public void testInterpolatedMedian5()
    {
        double val = Doubles.interpolatedMedian(Arrays.stream((new double[] { 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0 }))
                .boxed()
                .collect(Collectors.toList()));
        assertEquals(2.5, val, 1e-10);
    }
}
