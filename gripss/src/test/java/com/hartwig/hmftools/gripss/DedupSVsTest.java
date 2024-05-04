package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.IMPRECISE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;
import static com.hartwig.hmftools.gripss.GripssTestUtils.loadSvDataCache;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_TUMOR_AF;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.links.AlternatePath;
import com.hartwig.hmftools.gripss.links.Link;

import org.junit.Test;

public class DedupSVsTest
{
    private final GripssTestApp mGripss;
    private final SvDataCache mDataCache;
    private final FilterCache mFilterCache;
    private final DuplicateFinder mDuplicateFinder;

    public DedupSVsTest()
    {
        mGripss = new GripssTestApp();
        mDataCache = mGripss.DataCache;
        mFilterCache = mGripss.FilterCache;
        mDuplicateFinder = new DuplicateFinder(mDataCache, mFilterCache);
    }

    @Test
    public void testBasicDedup()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();

        commonAttributes.put(IMPRECISE, "true");

        SvData var1 = mGripss.createDel(CHR_1, 100, 1000, commonAttributes, commonAttributes);
        SvData var2 = mGripss.createDel(CHR_1, 100, 1000, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        AlternatePath altPath = new AlternatePath(var1.breakendStart(), var1.breakendEnd(), Lists.newArrayList(Link.from(var2)));

        // prioritisation is precise, then passing, then qual
        mDuplicateFinder.findDuplicateSVs(Lists.newArrayList(altPath));
        assertTrue(mDuplicateFinder.duplicateBreakends().contains(var1.breakendStart()));

        // test the reverse
        var1 = mGripss.createDel(CHR_1, 100, 1000, null, null);
        var2 = mGripss.createDel(CHR_1, 100, 1000, commonAttributes, commonAttributes);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        altPath = new AlternatePath(var1.breakendStart(), var1.breakendEnd(), Lists.newArrayList(Link.from(var2)));

        mDuplicateFinder.findDuplicateSVs(Lists.newArrayList(altPath));
        assertFalse(mDuplicateFinder.duplicateBreakends().contains(var1.breakendStart()));

        // passing vs not
        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(QUAL, DEFAULT_QUAL + 10); // higher than the other

        var1 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 100, 1000, POS_ORIENT, NEG_ORIENT, "",
                mGripss.GenotypeIds, null, null, tumorAttributes);

        var2 = mGripss.createDel(CHR_1, 100, 1000, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));
        altPath = new AlternatePath(var1.breakendStart(), var1.breakendEnd(), Lists.newArrayList(Link.from(var2)));

        mFilterCache.addBreakendFilter(var1.breakendStart(), MIN_TUMOR_AF);
        mDuplicateFinder.findDuplicateSVs(Lists.newArrayList(altPath));
        assertTrue(mDuplicateFinder.duplicateBreakends().contains(var1.breakendStart()));

        mDuplicateFinder.clear();

        mFilterCache.getBreakendFilters().remove(var1.breakendStart());
        mDuplicateFinder.findDuplicateSVs(Lists.newArrayList(altPath));
        assertFalse(mDuplicateFinder.duplicateBreakends().contains(var1.breakendStart()));
    }

    @Test
    public void testMultipleLinkDedup()
    {
        SvData var1 = mGripss.createDel(CHR_1, 100, 1000, null, null);
        SvData var2 = mGripss.createDel(CHR_1, 100, 400, null, null);
        SvData var3 = mGripss.createDel(CHR_1, 500, 1000, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        AlternatePath altPath = new AlternatePath(
                var1.breakendStart(), var1.breakendEnd(), Lists.newArrayList(Link.from(var2), Link.from(var3)));

        // prioritisation is precise, then passing, then qual
        mFilterCache.addBreakendFilter(var1.breakendStart(), MIN_TUMOR_AF);
        mDuplicateFinder.findDuplicateSVs(Lists.newArrayList(altPath));
        assertTrue(mDuplicateFinder.duplicateBreakends().contains(var1.breakendStart()));

        mDuplicateFinder.clear();

        // now with the main variant passing and one of the others not - still mark the variant as a dup and rescue to the others
        mFilterCache.getBreakendFilters().remove(var1.breakendStart());
        mFilterCache.addBreakendFilter(var3.breakendEnd(), MIN_TUMOR_AF);

        mDuplicateFinder.findDuplicateSVs(Lists.newArrayList(altPath));
        assertTrue(mDuplicateFinder.duplicateBreakends().contains(var1.breakendStart()));
        assertTrue(mDuplicateFinder.rescueBreakends().contains(var3.breakendEnd()));
    }
}