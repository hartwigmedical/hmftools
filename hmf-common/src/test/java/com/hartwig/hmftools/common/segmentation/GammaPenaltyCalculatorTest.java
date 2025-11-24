package com.hartwig.hmftools.common.segmentation;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GammaPenaltyCalculatorTest extends SegmentationTestBase
{
    @Test
    public void minimumCost()
    {
        GammaPenaltyCalculator gpc = new GammaPenaltyCalculator(28.0, true);
        assertEquals(0.01 * 28.0, gpc.getPenalty(d(1, 2, 3, 4)), 0.0001);
    }
}