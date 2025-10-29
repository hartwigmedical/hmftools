package com.hartwig.hmftools.cobalt.metrics;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class TargetRegionDataTest
{
    @Test
    public void constructorTest()
    {
        TargetRegionData data = new TargetRegionData("21", 2_000_001, 2_010_000);
        assertEquals(10_000, data.baseLength());
    }
}
