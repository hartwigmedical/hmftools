package com.hartwig.hmftools.common.sigs;

import static com.hartwig.hmftools.common.stats.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.stats.Percentiles.calcPercentileValues;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;

import static junit.framework.TestCase.assertEquals;


import org.junit.Test;

public class PercentilesTest
{
    @Test
    public void testPercentileCreation()
    {
        int slots = 11;
        double[] percentileValues = new double[slots];

        double[] values = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        calcPercentileValues(values, percentileValues);

        assertEquals(1.0, percentileValues[0], 0.001);
        assertEquals(1.4, percentileValues[1], 0.001);
        assertEquals(1.8, percentileValues[2], 0.001);
        assertEquals(2.2, percentileValues[3], 0.001);
        assertEquals(4.6, percentileValues[9], 0.001);
        assertEquals(5.0, percentileValues[10], 0.001);

        values = new double[51];

        for(int i = 0; i < values.length; ++i)
        {
            values[i] = i;
        }

        // will result in 5 values per slot
        calcPercentileValues(values, percentileValues);

        assertEquals(0, percentileValues[0], 0.001);
        assertEquals(5, percentileValues[1], 0.001);
        assertEquals(45, percentileValues[9], 0.001);
        assertEquals(50, percentileValues[10], 0.001);
    }

    @Test
    public void testPercentileSelection()
    {
        // first test with less values than slots
        double[] values = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double[] percentileValues = buildPercentiles(values);

        assertEquals(0, getPercentile(percentileValues, 0.5), 0.01);
        assertEquals(1.0, getPercentile(percentileValues, 6.0), 0.01);
        assertEquals(2.0, getPercentile(percentileValues, 12, true, 0), 0.01); // reports multiple of highest value
        assertEquals(0.4, getPercentile(percentileValues, 3.0), 0.01);
        assertEquals(0.7, getPercentile(percentileValues, 4.5), 0.01);

        // test with zeros for first X percentiles
        values = new double[] { 0.0, 0.0, 0.0, 4.0, 5.0, 6.0 };
        percentileValues = buildPercentiles(values);

        assertEquals(0.2, getPercentile(percentileValues, 0), 0.01);

        // test with zeros for all percentiles
        values = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        percentileValues = buildPercentiles(values);

        assertEquals(10.0, getPercentile(percentileValues, 1, true, 10.0), 0.01);

        // same value higher up
        values = new double[] { 1, 2, 3, 4, 4, 4, 4, 4, 5, 6, 7 };
        percentileValues = buildPercentiles(values);

        assertEquals(0.5, getPercentile(percentileValues, 4), 0.01);

        values = new double[] { 1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7 };
        percentileValues = buildPercentiles(values);

        assertEquals(0.8, getPercentile(percentileValues, 7), 0.01);
    }

}
