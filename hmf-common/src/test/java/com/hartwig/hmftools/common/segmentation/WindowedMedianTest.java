package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public class WindowedMedianTest
{
    @Test
    public void compareWithR()
    {
        check(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, 1);
        check(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, 3);
        check(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, 5);

        check(new double[] { 1.0, 3.0, 5.0, 2.0, 8.0, 7.0 }, new double[] { 1.0, 3.0, 5.0, 2.0, 8.0, 7.0 }, 1);
        check(new double[] { 1.0, 3.0, 5.0, 2.0, 8.0, 7.0 }, new double[] { 1.0, 3.0, 3.0, 5.0, 7.0, 7.0 }, 3);
        check(new double[] { 1.0, 3.0, 5.0, 2.0, 8.0, 7.0 }, new double[] { 1.0, 3.0, 3.0, 5.0, 8.0, 7.0 }, 5);

        check(new double[] { 1.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 }, new double[] { 1.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                9.0 }, 1);
        check(new double[] { 1.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 }, new double[] { 1.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                9.0 }, 3);
        check(new double[] { 1.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 }, new double[] { 1.0, 0.0, 1.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                9.0 }, 5);

        check(new double[] { 1.0, 0.0, -1.0, 2.0, 5.0, 6.0, 1.0, 8.0, 9.0 }, new double[] { 1.0, 0.0, 0.0, 2.0, 5.0, 5.0, 6.0, 8.0,
                9.0 }, 3);
        check(new double[] { 1.0, 0.0, -1.0, 2.0, 5.0, 6.0, 1.0, 8.0, 9.0 }, new double[] { 1.0, 0.0, 1.0, 2.0, 2.0, 5.0, 6.0, 8.0,
                9.0 }, 5);
    }

    private void check(double[] data, double[] expected, int windowSize)
    {
        double[] actual = new WindowedMedian(data, windowSize).getMedians();
        assertArrayEquals(expected, actual, 1e-10);
    }
}
