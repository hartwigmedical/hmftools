package com.hartwig.hmftools.cobalt.ratio;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;

import org.junit.Test;

public class DiploidRatioNormalizationTest
{
    private static final double EPSILON = 1e-10;

    @Test
    public void testCloseToZero()
    {
        final List<Double> input = Arrays.asList(0.0, 0.0, 0.002, 0.0, 0.0);

        final List<Double> output = new DiploidRatioNormalization(1.0, 5, 5, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1);
        assertRatio(input.get(1), output.get(1), 1);
        assertRatio(input.get(2), output.get(2), 1);
        assertRatio(input.get(3), output.get(3), 1);
        assertRatio(input.get(4), output.get(4), 1);
    }

    @Test
    public void testMaxWindowDistance()
    {
        final List<Double> input = Arrays.asList(1.0, 1.5, -1.0, 1.1, 1.2);

        final List<Double> output = new DiploidRatioNormalization(1.0, 2, 1, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1.25);
        assertRatio(input.get(1), output.get(1), 1.1);
        assertRatio(input.get(2), output.get(2), 1.0);
        assertRatio(input.get(3), output.get(3), 1.2);
        assertRatio(input.get(4), output.get(4), 1.15);
    }

    @Test
    public void testMinCoverage()
    {
        final List<Double> input = Arrays.asList(1.0, 1.5, 2.0, -1.0, -1.0);

        final List<Double> output = new DiploidRatioNormalization(1.0, 1, 3, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1.0);
        assertRatio(input.get(1), output.get(1), 1.5);
        assertRatio(input.get(2), output.get(2), 1.0);
        assertRatio(input.get(3), output.get(3), 1.0);
        assertRatio(input.get(4), output.get(4), 1.0);
    }

    private static void assertRatio(final double input, final double output, double median)
    {
        assertEquals(input / median, output, EPSILON);
    }
}
