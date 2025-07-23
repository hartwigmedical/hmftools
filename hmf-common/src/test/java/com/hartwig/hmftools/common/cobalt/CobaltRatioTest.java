package com.hartwig.hmftools.common.cobalt;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CobaltRatioTest
{
    @Test
    public void realignTest()
    {
        CobaltRatio ratio = new CobaltRatio("chr1", 1000, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8);
        int newPosition = 2000;
        CobaltRatio realigned = ratio.realign(newPosition);
        assertEquals(newPosition, realigned.position());
        assertEquals(ratio.chromosome(), realigned.chromosome());
        assertEquals(ratio.referenceReadDepth(), realigned.referenceReadDepth(),0.00001);
        assertEquals(ratio.referenceGCRatio(), realigned.referenceGCRatio(), 0.0001);
        assertEquals(ratio.referenceGcContent(), realigned.referenceGcContent(), 0.0001);
        assertEquals(ratio.referenceGCDiploidRatio(), realigned.referenceGCDiploidRatio(), 0.0001);
        assertEquals(ratio.tumorReadDepth(), realigned.tumorReadDepth(),0.00001);
        assertEquals(ratio.tumorGcContent(), realigned.tumorGcContent(),0.00001);
        assertEquals(ratio.tumorGCRatio(), realigned.tumorGCRatio(),0.00001);
    }
}
