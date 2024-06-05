package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BAQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.IMPRECISE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSgl;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;
import static com.hartwig.hmftools.gripss.GripssTestUtils.loadSvDataCache;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.links.Link;
import com.hartwig.hmftools.gripss.links.LinkStore;

import org.junit.Test;

public class DedupSinglesTest
{
    private final GripssTestApp mGripss;
    private final SvDataCache mDataCache;
    private final FilterCache mFilterCache;
    private final DuplicateFinder mDuplicateFinder;

    private static final Map<String, Object> DEFAULT_ATTRIBUTES = Maps.newHashMap();

    public DedupSinglesTest()
    {
        mGripss = new GripssTestApp();
        mDataCache = mGripss.DataCache;
        mFilterCache = mGripss.FilterCache;
        mDuplicateFinder = new DuplicateFinder(mDataCache, mFilterCache);

        DEFAULT_ATTRIBUTES.put(CIPOS, new int[] {-10,10});
    }

    @Test
    public void testSingleDedup()
    {
        SvData var1 = mGripss.createDel(CHR_1, 101, 1000, DEFAULT_ATTRIBUTES,  DEFAULT_ATTRIBUTES);

        SvData var2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(var2.breakendStart()));
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(var1.breakendStart()));
    }

    @Test
    public void testKeepSingleOverPairBecauseOfQual()
    {
        SvData var1 = mGripss.createDel(CHR_1, 101, 1000, DEFAULT_ATTRIBUTES,  DEFAULT_ATTRIBUTES);

        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(GRIDSS_BAQ, DEFAULT_QUAL + 10);

        SvData var2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, tumorAttributes);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(var1.breakendStart()));
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(var2.breakendStart()));
    }

    @Test
    public void testKeepPairBecauseOfLink()
    {
        SvData var1 = mGripss.createDel(CHR_1, 101, 1000, DEFAULT_ATTRIBUTES,  DEFAULT_ATTRIBUTES);

        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(GRIDSS_BAQ, DEFAULT_QUAL + 10);

        SvData var2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, tumorAttributes);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        LinkStore linkStore = new LinkStore();
        linkStore.addLink(var1.breakendStart(), Link.from(var1));

        mDuplicateFinder.findDuplicateSingles(linkStore);
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(var1.breakendStart()));
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(var2.breakendStart()));
    }

    @Test
    public void testKeepSingleOverPairBecauseOfPairFailing()
    {
        SvData var1 = mGripss.createDel(CHR_1, 101, 1000, DEFAULT_ATTRIBUTES,  DEFAULT_ATTRIBUTES);

        SvData var2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        mFilterCache.addBreakendFilter(var1.breakendStart(), FilterType.MIN_TUMOR_AF);

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(var1.breakendStart()));
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(var2.breakendStart()));
    }

    @Test
    public void testOutsideRange()
    {
        SvData var1 = mGripss.createDel(CHR_1, 110, 1000, DEFAULT_ATTRIBUTES,  DEFAULT_ATTRIBUTES);

        SvData var2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 89, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, var2));

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertTrue(mDuplicateFinder.duplicateSglBreakends().isEmpty());
    }

    @Test
    public void testMultipleSinglesAndPair()
    {
        SvData var1 = mGripss.createDel(CHR_1, 103, 1000, DEFAULT_ATTRIBUTES,  DEFAULT_ATTRIBUTES);

        SvData sgl1 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        SvData sgl2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 102, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, sgl1, sgl2));

        mFilterCache.addBreakendFilter(var1.breakendStart(), FilterType.MIN_TUMOR_AF);

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(var1.breakendStart()));
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(sgl1.breakendStart()));
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(sgl2.breakendStart()));
    }

    @Test
    public void testChoosePassingSingle()
    {
        SvData sgl1 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        SvData sgl2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 102, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(sgl1, sgl2));

        mFilterCache.addBreakendFilter(sgl2.breakendStart(), FilterType.MIN_TUMOR_AF);

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(sgl1.breakendStart()));
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(sgl2.breakendStart()));
    }

    @Test
    public void testChooseHighestQuality()
    {
        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(GRIDSS_BAQ, 10);

        SvData sgl1 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, tumorAttributes);

        SvData sgl2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 102, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(sgl1, sgl2));

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(sgl1.breakendStart()));
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(sgl2.breakendStart()));
    }

    @Test
    public void testDifferentOrientation()
    {
        SvData sgl1 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 100, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        SvData sgl2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 102, NEG_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(sgl1, sgl2));

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertTrue(mDuplicateFinder.duplicateSglBreakends().isEmpty());
    }

    @Test
    public void testDoNotUseConfidenceIntervalOfImprecise()
    {
        SvData sgl1 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 90, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        Map<String, Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(IMPRECISE, "true");

        SvData sgl2 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 101, NEG_ORIENT, "",
                mGripss.GenotypeIds, commonAttributes, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(sgl1, sgl2));

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertTrue(mDuplicateFinder.duplicateSglBreakends().isEmpty());
    }

    @Test
    public void testCanStillMatchAgainstImprecise()
    {
        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(QUAL, 6000);

        Map<String, Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(CIPOS, new int[] {-10, 10});
        commonAttributes.put(IMPRECISE, "true");

        SvData var1 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 100, 2000, POS_ORIENT, NEG_ORIENT, "",
                mGripss.GenotypeIds, commonAttributes, null, tumorAttributes);

        SvData sgl1 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 90, POS_ORIENT, "",
                mGripss.GenotypeIds, DEFAULT_ATTRIBUTES, null, null);

        loadSvDataCache(mDataCache, Lists.newArrayList(var1, sgl1));

        mDuplicateFinder.findDuplicateSingles(new LinkStore());
        assertFalse(mDuplicateFinder.duplicateSglBreakends().contains(var1.breakendStart()));
        assertTrue(mDuplicateFinder.duplicateSglBreakends().contains(sgl1.breakendStart()));
    }
}
