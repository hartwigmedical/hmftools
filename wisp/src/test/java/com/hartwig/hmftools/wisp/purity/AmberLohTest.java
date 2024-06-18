package com.hartwig.hmftools.wisp.purity;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.wisp.purity.loh.LodCalcs;

import org.junit.Test;

public class AmberLohTest
{
    @Test
    public void testAmberLohLod()
    {
        double lod = LodCalcs.calcLimitOfDetection(0.01, 2, 50);
        assertEquals(0.45, lod, 0.01);

        lod = LodCalcs.calcLimitOfDetection(0.01, 2, 100);
        assertEquals(0.32, lod, 0.01);

        lod = LodCalcs.calcLimitOfDetection(0.01, 2, 1000);
        assertEquals(0.10, lod, 0.01);

        lod = LodCalcs.calcLimitOfDetection(0.01, 2, 10000);
        assertEquals(0.03, lod, 0.01);

        lod = LodCalcs.calcLimitOfDetection(0.011, 2, 100000);
        assertEquals(0.01, lod, 0.001);

        lod = LodCalcs.calcLimitOfDetection(0.0033, 2, 1000000);
        assertEquals(0.0038, lod, 0.001);
    }
}
