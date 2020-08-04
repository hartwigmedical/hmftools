package com.hartwig.hmftools.common.sigs;

import static com.hartwig.hmftools.common.sigs.Percentiles.calcPercentileValues;

import static junit.framework.TestCase.assertEquals;


import org.junit.Test;

public class PercentilesTest
{
    @Test
    public void testPercentileSplits()
    {
        int slots = 11;
        double[] percentileValues = new double[slots];

        double[] values = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        calcPercentileValues(values, percentileValues);

        assertEquals(1.0, percentileValues[0], 0.001);
        assertEquals(1.0, percentileValues[1], 0.001);
        assertEquals(1.8, percentileValues[2], 0.001);
        assertEquals(2.0, percentileValues[3], 0.001);
        assertEquals(5.0, percentileValues[9], 0.001);
        assertEquals(5.0, percentileValues[10], 0.001);

        values = new double[50];

        for(int i = 0; i < values.length; ++i)
        {
            values[i] = i;
        }

        calcPercentileValues(values, percentileValues);

        assertEquals(1.8, percentileValues[0], 0.001);
        assertEquals(6.3, percentileValues[1], 0.001);
        assertEquals(42.7, percentileValues[9], 0.001);
        assertEquals(47.2, percentileValues[10], 0.001);
    }
}
