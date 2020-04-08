package com.hartwig.hmftools.isofox;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;

import org.junit.Test;

public class GcBiasTest
{
    @Test
    public void testGcRatioCalcs()
    {
        GcRatioCounts gcCounts = new GcRatioCounts();

        final int[] gcRatioIndex = {-1, -1};
        final double[] counts = {1, 0};

        double gcRatio = 0;
        gcCounts.determineRatioData(gcRatio, gcRatioIndex, counts);
        assertEquals(0, gcRatioIndex[0]);
        assertEquals(1, gcRatioIndex[1]);
        assertEquals(1, counts[0], 0.001);
        assertEquals(0, counts[1], 0.001);

        gcRatio = 0.505;
        gcCounts.determineRatioData(gcRatio, gcRatioIndex, counts);
        assertEquals(50, gcRatioIndex[0]);
        assertEquals(51, gcRatioIndex[1]);
        assertEquals(0.5, counts[0], 0.001);
        assertEquals(0.5, counts[1], 0.001);

        gcRatio = 0.9925;
        gcCounts.determineRatioData(gcRatio, gcRatioIndex, counts);
        assertEquals(99, gcRatioIndex[0]);
        assertEquals(100, gcRatioIndex[1]);
        assertEquals(0.75, counts[0], 0.001);
        assertEquals(0.25, counts[1], 0.001);

        gcRatio = 1.0;
        gcCounts.determineRatioData(gcRatio, gcRatioIndex, counts);
        assertEquals(100, gcRatioIndex[0]);
        assertEquals(1, counts[0], 0.001);
    }


}
