package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

public class RunmedTest
{
    @Test
    public void compareWithR()
    {
        check(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, 1);
        check(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, 3);
        check(new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 }, 5);

        check(new double[] { 1.0, 1.0, 2.0, 2.0, 3.0, 3.0 }, new double[] { 1.0, 1.0, 2.0, 2.0, 3.0, 3.0 }, 3);

        check(new double[] { 1.0, 3.0, 5.0, 2.0, 8.0, 7.0, 8.0, 3.0, 12.0 }, new double[] { 3.0, 3.0, 3.0, 5.0, 7.0, 8.0, 7.0, 8.0,
                10.0 }, 3);
        check(new double[] { 1.0, 3.0, 5.0, 2.0, 8.0, 7.0, 8.0, 3.0, 12.0 }, new double[] { 3.0, 3.0, 3.0, 5.0, 7.0, 7.0, 8.0, 8.0,
                8.0 }, 5);

        check(new double[] { 3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0, 9.0 }, new double[] { 3.0, 3.0, 1.0, 4.0, 5.0, 5.0, 6.0,
                5.0, 5.0, 5.0, 5.0, 5.0 }, 3);
        check(new double[] { 3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0, 9.0 }, new double[] { 3.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0,
                5.0, 5.0, 5.0, 5.0, 5.0 }, 5);
        check(new double[] { 3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0, 9.0 }, new double[] { 3.0, 3.0, 3.0, 3.0, 4.0, 5.0, 5.0,
                5.0, 5.0, 5.0, 5.0, 5.0 }, 7);

        // Negative values
        check(new double[] { 1.0, 1.0, -5.0, -6.0, 1.0, 2.0 }, new double[] { 1.0, 1.0, -5.0, -6.0, 1.0, 2.0 }, 1);
        check(new double[] { 1.0, 1.0, -5.0, -6.0, 1.0, 2.0 }, new double[] { 1.0, 1.0, -5.0, -5.0, 1.0, 2.0 }, 3);
        check(new double[] { 1.0, 1.0, -5.0, -6.0, 1.0, 2.0 }, new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, 5);
    }

    @Test
    public void empty()
    {
        check(new double[] {}, new double[] {}, 1);
    }

    @Test
    public void singleElement()
    {
        check(new double[] { 42.0 }, new double[] { 42.0 }, 1);
    }

    @Test
    public void testRunmedInvalidWindowSize()
    {
        double[] data = new double[] { 1.0, 2.0, 3.0 };

        try
        {
            new Runmed().runmed(data, 2);
            fail("Expected IllegalArgumentException");
        }
        catch(IllegalArgumentException e)
        {
            // Expected
        }

        try
        {
            new Runmed().runmed(data, 0);
            fail("Expected IllegalArgumentException");
        }
        catch(IllegalArgumentException e)
        {
            // Expected
        }

        try
        {
            new Runmed().runmed(data, -1);
            fail("Expected IllegalArgumentException");
        }
        catch(IllegalArgumentException e)
        {
            // Expected
        }
    }

    private void check(double[] data, double[] expected, int windowSize)
    {
        double[] actual = new Runmed().runmed(data, windowSize);
        assertArrayEquals(expected, actual, 1e-10);
    }
}