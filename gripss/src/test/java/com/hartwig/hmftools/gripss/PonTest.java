package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_CIPOS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IHOMPOS;
import static com.hartwig.hmftools.gripss.pon.PonCombiner.mergeSglRegions;
import static com.hartwig.hmftools.gripss.pon.PonCombiner.mergeSvRegions;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.pon.PonCache;
import com.hartwig.hmftools.gripss.pon.PonSglRegion;
import com.hartwig.hmftools.gripss.pon.PonSvRegion;

import org.junit.Test;

public class PonTest
{
    private final GripssTestApp mGripss;
    private final PonCache mPonCache;

    public PonTest()
    {
        mGripss = new GripssTestApp();
        mPonCache = new PonCache(2, null, null, false);
    }

    @Test
    public void testSvPonMatches()
    {
        int startPos = 100;
        int endPos = 1000;
        int ponBuffer = 10;
        int ponGap = 1000;

        for(int i = 0; i < 16; ++i)
        {
            mPonCache.addPonSvRegion(
                    CHR_1, new BaseRegion(startPos, startPos + ponBuffer), POS_ORIENT,
                    new ChrBaseRegion(CHR_1, endPos, endPos + ponBuffer), NEG_ORIENT, 1);

            startPos += ponGap;
            endPos += ponGap;
        }

        SvData var1 = mGripss.createDel(CHR_1, 100, 500, null, null);
        assertFalse(mPonCache.getPonCount(var1) > 0);

        var1 = mGripss.createDel(CHR_1, 100, 1005, null, null);
        assertTrue(mPonCache.getPonCount(var1) > 0);

        // must match on orientation
        var1 = mGripss.createDup(CHR_1, 2100, 3005, null, null);
        assertFalse(mPonCache.getPonCount(var1) > 0);

        // requires margin to match
        var1 = mGripss.createDel(CHR_1, 3098, 4012, null, null);
        assertTrue(mPonCache.getPonCount(var1) > 0);

        Map<String,Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(VT_IHOMPOS, new int[] {-10, 10});

        // match within inexact homology
        var1 = mGripss.createDel(CHR_1, 4088, 5022, commonAttributes, commonAttributes);
        assertTrue(mPonCache.getPonCount(var1) > 0);

        var1 = mGripss.createDel(CHR_1, 5087, 6023, commonAttributes, commonAttributes);
        assertFalse(mPonCache.getPonCount(var1) > 0);

        // now test with ranges that are smaller and larger and inexact variations
        mPonCache.clear();

        mPonCache.addPonSvRegion(
                CHR_1, new BaseRegion(100, 102), POS_ORIENT,
                new ChrBaseRegion(CHR_1, 200, 210), NEG_ORIENT, 1);

        mPonCache.addPonSvRegion(
                CHR_1, new BaseRegion(101, 102), POS_ORIENT,
                new ChrBaseRegion(CHR_1, 150, 151), NEG_ORIENT, 1);

        mPonCache.addPonSvRegion(
                CHR_1, new BaseRegion(102, 104), POS_ORIENT,
                new ChrBaseRegion(CHR_1, 200, 210), NEG_ORIENT, 1);

        // first SV will not match and move the index ahead of where the second one needs to match
        var1 = mGripss.createDel(CHR_1, 108, 305, null, null);
        assertFalse(mPonCache.getPonCount(var1) > 0);

        commonAttributes.put(VT_CIPOS, new int[] {-15, 15});
        commonAttributes.put(VT_IHOMPOS, new int[] {-20, 20});

        var1 = mGripss.createDel(CHR_1, 110, 205, commonAttributes, commonAttributes);
        assertTrue(mPonCache.getPonCount(var1) > 0);

    }

    @Test
    public void testSglPonMatches()
    {
        int startPos = 100;
        int ponBuffer = 10;
        int ponGap = 1000;

        for(int i = 0; i < 16; ++i)
        {
            mPonCache.addPonSglRegion(CHR_1, new BaseRegion(startPos, startPos + ponBuffer), POS_ORIENT, 1);

            startPos += ponGap;
        }

        SvData var = GripssTestUtils.createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, null, null, null);

        assertTrue(mPonCache.getPonCount(var) > 0);

        // requires margin to match
        var = GripssTestUtils.createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 1098, POS_ORIENT, "",
                mGripss.GenotypeIds, null, null, null);

        assertTrue(mPonCache.getPonCount(var) > 0);

        var = GripssTestUtils.createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 2113, POS_ORIENT, "",
                mGripss.GenotypeIds, null, null, null);

        assertFalse(mPonCache.getPonCount(var) > 0);

        Map<String,Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(VT_IHOMPOS, new int[] {-10, 0});

        // match within inexact homology
        var = GripssTestUtils.createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 3088, POS_ORIENT, "",
                mGripss.GenotypeIds, commonAttributes, null, null);

        assertTrue(mPonCache.getPonCount(var) > 0);

        var = GripssTestUtils.createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 8122, POS_ORIENT, "",
                mGripss.GenotypeIds, commonAttributes, null, null);

        assertTrue(mPonCache.getPonCount(var) > 0);
    }

    @Test
    public void testPonSvRegionMerge()
    {
        List<PonSvRegion> combinedRegions = Lists.newArrayList();

        combinedRegions.add(new PonSvRegion(
                new BaseRegion(100, 105), POS_ORIENT,
                new ChrBaseRegion(CHR_2, 100, 110), NEG_ORIENT, 10));

        // will be combine with 2
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(98, 103), POS_ORIENT,
                new ChrBaseRegion(CHR_2, 102, 104), NEG_ORIENT, 10));

        // not combined, diff orientation
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(102, 107), NEG_ORIENT,
                new ChrBaseRegion(CHR_2, 100, 110), NEG_ORIENT, 10));

        // not combined, latter breakend doesn't overlap
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(104, 109), POS_ORIENT,
                new ChrBaseRegion(CHR_2, 120, 125), NEG_ORIENT, 10));

        // latter start but still overlaps, will be combined
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(105, 106), POS_ORIENT,
                new ChrBaseRegion(CHR_2, 105, 115), NEG_ORIENT, 10));

        // new region
        combinedRegions.add(new PonSvRegion(
                new BaseRegion(200, 205), POS_ORIENT,
                new ChrBaseRegion(CHR_2, 100, 105), NEG_ORIENT, 10));

        combinedRegions.add(new PonSvRegion(
                new BaseRegion(202, 207), POS_ORIENT,
                new ChrBaseRegion(CHR_2, 95, 100), NEG_ORIENT, 10));

        mergeSvRegions(CHR_1, combinedRegions);

        assertEquals(4, combinedRegions.size());

        // check combined regions
        assertEquals(98, combinedRegions.get(0).RegionStart.start());
        assertEquals(106, combinedRegions.get(0).RegionStart.end());
        assertEquals(100, combinedRegions.get(0).RegionEnd.start());
        assertEquals(115, combinedRegions.get(0).RegionEnd.end());

        // unchanged
        assertEquals(102, combinedRegions.get(1).RegionStart.start());
        assertEquals(107, combinedRegions.get(1).RegionStart.end());
        assertEquals(100, combinedRegions.get(1).RegionEnd.start());
        assertEquals(110, combinedRegions.get(1).RegionEnd.end());

        // unchanged
        assertEquals(104, combinedRegions.get(2).RegionStart.start());
        assertEquals(109, combinedRegions.get(2).RegionStart.end());
        assertEquals(120, combinedRegions.get(2).RegionEnd.start());
        assertEquals(125, combinedRegions.get(2).RegionEnd.end());

        // merged
        assertEquals(200, combinedRegions.get(3).RegionStart.start());
        assertEquals(207, combinedRegions.get(3).RegionStart.end());
        assertEquals(95, combinedRegions.get(3).RegionEnd.start());
        assertEquals(105, combinedRegions.get(3).RegionEnd.end());
    }

    @Test
    public void testPonSglRegionMerge()
    {
        List<PonSglRegion> combinedRegions = Lists.newArrayList();

        combinedRegions.add(new PonSglRegion(
                new BaseRegion(100, 105), POS_ORIENT, 10));

        // will be combine with 2
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(98, 103), POS_ORIENT, 10));

        // not combined, diff orientation
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(100, 105), NEG_ORIENT, 10));

        // latter start but still overlaps, will be combined
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(105, 106), POS_ORIENT, 10));

        // new region
        combinedRegions.add(new PonSglRegion(
                new BaseRegion(200, 205), NEG_ORIENT, 10));

        combinedRegions.add(new PonSglRegion(
                new BaseRegion(203, 210), NEG_ORIENT, 10));

        mergeSglRegions(CHR_1, combinedRegions);

        assertEquals(3, combinedRegions.size());

        // check combined regions
        assertEquals(98, combinedRegions.get(0).Region.start());
        assertEquals(106, combinedRegions.get(0).Region.end());

        // unchanged
        assertEquals(100, combinedRegions.get(1).Region.start());
        assertEquals(105, combinedRegions.get(1).Region.end());

        // merged
        assertEquals(200, combinedRegions.get(2).Region.start());
        assertEquals(210, combinedRegions.get(2).Region.end());
    }
}
