package com.hartwig.hmftools.common.cobalt;

import static java.util.List.of;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Assert;
import org.junit.Test;

public class CobaltRatioTest
{
    CobaltRatio ratio = new CobaltRatio("chr1", 1001, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8);

    @Test
    public void realignTest()
    {
        int newPosition = 2001;
        CobaltRatio realigned = ratio.realign(newPosition);
        assertEquals(newPosition, realigned.position());
        assertEquals(ratio.chromosome(), realigned.chromosome());
        assertEquals(ratio.referenceReadDepth(), realigned.referenceReadDepth(), 0.00001);
        assertEquals(ratio.referenceGCRatio(), realigned.referenceGCRatio(), 0.0001);
        assertEquals(ratio.referenceGcContent(), realigned.referenceGcContent(), 0.0001);
        assertEquals(ratio.referenceGCDiploidRatio(), realigned.referenceGCDiploidRatio(), 0.0001);
        assertEquals(ratio.tumorReadDepth(), realigned.tumorReadDepth(), 0.00001);
        assertEquals(ratio.tumorGcContent(), realigned.tumorGcContent(), 0.00001);
        assertEquals(ratio.tumorGCRatio(), realigned.tumorGCRatio(), 0.00001);
    }

    @Test
    public void maskTest()
    {
        CobaltRatio masked = ratio.mask();
        assertEquals(ratio.position(), masked.position());
        assertEquals(ratio.chromosome(), masked.chromosome());
        assertEquals(ratio.referenceReadDepth(), masked.referenceReadDepth(), 0.00001);
        assertEquals(-1.0, masked.referenceGCRatio(), 0.0001);
        assertEquals(-1.0, masked.referenceGcContent(), 0.0001);
        assertEquals(-1.0, masked.referenceGCDiploidRatio(), 0.0001);
        assertEquals(ratio.tumorReadDepth(), masked.tumorReadDepth(), 0.00001);
        assertEquals(-1.0, masked.tumorGcContent(), 0.00001);
        assertEquals(-1.0, masked.tumorGCRatio(), 0.00001);
    }

    @Test
    public void windowTest()
    {
        assertEquals(new ChrBaseRegion(ratio.chromosome(), 1001, 2000), ratio.window());
    }

    @Test
    public void overlapsEmptyTest()
    {
        Assert.assertTrue(ratio.findWindowOverlaps(Collections.emptyList()).isEmpty());
    }

    @Test
    public void noOverlapsTest()
    {
        List<ChrBaseRegion> regions = of(gr(100, 200), gr(400, 500), gr(2100, 2300));
        Assert.assertTrue(ratio.findWindowOverlaps(regions).isEmpty());
    }

    @Test
    public void overlapMissByOne()
    {
        int oneBefore = ratio.position() - 1;
        int oneAfter = ratio.window().end() + 1;
        List<ChrBaseRegion> regions = of(gr(oneBefore - 100, oneBefore), gr(oneAfter, oneAfter + 100));
        Assert.assertTrue(ratio.findWindowOverlaps(regions).isEmpty());
    }

    @Test
    public void overlapWithFirstBase()
    {
        ChrBaseRegion region = gr(ratio.position() - 100, ratio.position());
        List<ChrBaseRegion> overlaps = ratio.findWindowOverlaps(of(region));
        assertEquals(1, overlaps.size());
        assertEquals(region, overlaps.get(0));
    }

    @Test
    public void overlapsWithLastBase()
    {
        ChrBaseRegion region = gr(ratio.window().end(), ratio.window().end() + 100);
        List<ChrBaseRegion> overlaps = ratio.findWindowOverlaps(of(region));
        assertEquals(1, overlaps.size());
        assertEquals(region, overlaps.get(0));
    }

    @Test
    public void multipleOverlaps()
    {
        ChrBaseRegion region0 = gr(850, 950);
        ChrBaseRegion region1 = gr(950, 1050);
        ChrBaseRegion region2 = gr(1150, 1250);
        ChrBaseRegion region3 = gr(1750, 1850);
        ChrBaseRegion region4 = gr(1950, 2050);
        ChrBaseRegion region5 = gr(2150, 2250);
        List<ChrBaseRegion> overlaps = ratio.findWindowOverlaps(of(region0, region1, region2, region3, region4, region5));
        assertEquals(of(region1, region2, region3, region4), overlaps);
    }

    private ChrBaseRegion gr(int start, int end)
    {
        return new ChrBaseRegion(ratio.chromosome(), start, end);
    }
}