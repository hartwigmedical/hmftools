package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

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
        commonAttributes.put(IHOMPOS, new int[] {-10, 10});

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

        commonAttributes.put(CIPOS, new int[] {-15, 15});
        commonAttributes.put(IHOMPOS, new int[] {-20, 20});

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
        commonAttributes.put(IHOMPOS, new int[] {-10, 0});

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
}
