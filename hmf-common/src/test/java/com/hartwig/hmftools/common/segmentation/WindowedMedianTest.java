package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

public class WindowedMedianTest extends SegmentationTestBase
{
    @Test
    public void compareWithR()
    {
        check(d(1, 2, 3, 4, 5), d(1, 2, 3, 4, 5), 1);
        check(d(1, 2, 3, 4, 5), d(1, 2, 3, 4, 5), 3);
        check(d(1, 2, 3, 4, 5), d(1, 2, 3, 4, 5), 5);

        check(d(1, 3, 5, 2, 8, 7), d(1, 3, 5, 2, 8, 7), 1);
        check(d(1, 3, 5, 2, 8, 7), d(1, 3, 3, 5, 7, 7), 3);
        check(d(1, 3, 5, 2, 8, 7), d(1, 3, 3, 5, 8, 7), 5);

        check(d(1, 0, 0, 4, 5, 6, 7, 8, 9), d(1, 0, 0, 4, 5, 6, 7, 8, 9), 1);
        check(d(1, 0, 0, 4, 5, 6, 7, 8, 9), d(1, 0, 0, 4, 5, 6, 7, 8, 9), 3);
        check(d(1, 0, 0, 4, 5, 6, 7, 8, 9), d(1, 0, 1, 4, 5, 6, 7, 8, 9), 5);

        check(d(1, 0, -1, 2, 5, 6, 1, 8, 9), d(1, 0, 0, 2, 5, 5, 6, 8, 9), 3);
        check(d(1, 0, -1, 2, 5, 6, 1, 8, 9), d(1, 0, 1, 2, 2, 5, 6, 8, 9), 5);
    }

    private void check(double[] data, double[] expected, int windowSize)
    {
        double[] actual = new WindowedMedian(data, windowSize).getMedians();
        assertArrayEquals(expected, actual, 1e-10);
    }
}
