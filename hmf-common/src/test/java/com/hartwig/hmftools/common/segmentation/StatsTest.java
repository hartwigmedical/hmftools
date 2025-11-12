package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class StatsTest extends SegmentationTestBase
{
    private static final double DELTA = 0.001;

    @Test
    public void medianAbsoluteDeviationTest()
    {
        assertEquals(0, Stats.medianAbsoluteDeviation(d(1)), DELTA);
        assertEquals(2.9652, Stats.medianAbsoluteDeviation(d(1, 3, 5)), DELTA);
        assertEquals(3.7064, Stats.medianAbsoluteDeviation(d(1, 3, 5, 2, 8, 7)), DELTA);
    }
}