package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IHOMPOS;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.gripss.common.SvData;

import org.junit.Test;

public class PonTest
{
    private final GripssTestApp mGripss;
    private final PonCache mPonCache;

    public PonTest()
    {
        mGripss = new GripssTestApp();
        mPonCache = new PonCache(2, null, null);
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
                    new ChrBaseRegion(CHR_1, startPos, startPos + ponBuffer), POS_ORIENT,
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
}
