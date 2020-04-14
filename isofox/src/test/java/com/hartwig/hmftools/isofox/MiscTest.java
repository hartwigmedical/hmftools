package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.common.RnaUtils.calcPercentileValues;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.calcEffectiveLength;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class MiscTest
{
    @Test
    public void testEffectiveLength()
    {
        final List<int[]> fragmentLengthData = Lists.newArrayList();

        fragmentLengthData.add(new int[]{100, 1});
        fragmentLengthData.add(new int[]{200, 2});
        fragmentLengthData.add(new int[]{300, 1});

        assertEquals(800, calcEffectiveLength(1000, fragmentLengthData), 0.001);

        assertEquals(83.33, calcEffectiveLength(250, fragmentLengthData), 0.1);
    }

    @Test
    public void testPercentileSplits()
    {
        int slots = 11;
        double[] percentileValues = new double[slots];

        List<Double> values = Lists.newArrayList(1.0, 2.0, 3.0, 4.0, 5.0);
        calcPercentileValues(values, percentileValues);

        assertEquals(1.0, percentileValues[0], 0.001);
        assertEquals(1.0, percentileValues[1], 0.001);
        assertEquals(1.8, percentileValues[2], 0.001);
        assertEquals(2.0, percentileValues[3], 0.001);
        assertEquals(5.0, percentileValues[9], 0.001);
        assertEquals(5.0, percentileValues[10], 0.001);

        values.clear();

        for(int i = 0; i < 50; ++i)
        {
            values.add((double)i);
        }

        calcPercentileValues(values, percentileValues);

        assertEquals(1.8, percentileValues[0], 0.001);
        assertEquals(6.3, percentileValues[1], 0.001);
        assertEquals(42.7, percentileValues[9], 0.001);
        assertEquals(47.2, percentileValues[10], 0.001);

    }
}
