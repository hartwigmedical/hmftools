package com.hartwig.hmftools.common.cobalt;

import static java.util.List.of;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

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
    public void windowTest()
    {
        assertEquals(new BaseRegion(1001, 2000), ratio.window());
    }

    @Test
    public void overlapsEmptyTest()
    {
        Assert.assertTrue(ratio.findWindowOverlaps(Collections.emptyList()).isEmpty());
    }

    @Test
    public void noOverlapsTest()
    {
        List<BaseRegion> regions = of(br(100, 200), br(400, 500), br(2100, 2300));
        Assert.assertTrue(ratio.findWindowOverlaps(regions).isEmpty());
    }

    @Test
    public void overlapMissByOne()
    {
        int oneBefore = ratio.position() - 1;
        int oneAfter = ratio.window().end() + 1;
        List<BaseRegion> regions = of(br(oneBefore - 100, oneBefore), br(oneAfter, oneAfter + 100));
        Assert.assertTrue(ratio.findWindowOverlaps(regions).isEmpty());
    }

    @Test
    public void overlapWithFirstBase()
    {
        BaseRegion region = br(ratio.position() - 100, ratio.position());
        List<BaseRegion> overlaps = ratio.findWindowOverlaps(of(region));
        assertEquals(1, overlaps.size());
        assertEquals(region, overlaps.get(0));
    }

    @Test
    public void overlapsWithLastBase()
    {
        BaseRegion region = br(ratio.window().end(), ratio.window().end() + 100);
        List<BaseRegion> overlaps = ratio.findWindowOverlaps(of(region));
        assertEquals(1, overlaps.size());
        assertEquals(region, overlaps.get(0));
    }

    @Test
    public void multipleOverlaps()
    {
        BaseRegion region0 = br(850, 950);
        BaseRegion region1 = br(950, 1050);
        BaseRegion region2 = br(1150, 1250);
        BaseRegion region3 = br(1750, 1850);
        BaseRegion region4 = br(1950, 2050);
        BaseRegion region5 = br(2150, 2250);
        List<BaseRegion> overlaps = ratio.findWindowOverlaps(of(region0, region1, region2, region3, region4, region5));
        assertEquals(of(region1, region2, region3, region4), overlaps);
    }

    private BaseRegion br(int start, int end)
    {
        return new BaseRegion(start, end);
    }
}