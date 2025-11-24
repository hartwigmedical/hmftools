package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class FixedPenaltyTest extends SegmentationTestBase
{
    @Test
    public void minimumCost()
    {
        FixedPenalty gpc = new FixedPenalty(28.0);
        assertEquals(28.0, gpc.getPenalty(d(1, 2, 3, 4)), 0.0001);
    }
}