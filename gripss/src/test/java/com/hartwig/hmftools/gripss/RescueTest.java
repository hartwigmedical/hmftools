package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.LINE_INSERT_SEQ_A;
import static com.hartwig.hmftools.gripss.GripssTestUtils.LINE_INSERT_SEQ_T;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_TUMOR_AF;

import static junit.framework.TestCase.assertTrue;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.links.LinkRescue;
import com.hartwig.hmftools.gripss.links.LinkStore;

import org.junit.Test;

public class RescueTest
{
    private final GripssTestApp mGripss;
    private final LinkRescue mLinkRescue;

    public RescueTest()
    {
        mGripss = new GripssTestApp();
        mLinkRescue = new LinkRescue();
    }

    private void setFiltered(final Breakend breakend) { mGripss.FilterCache.addBreakendFilter(breakend, MIN_TUMOR_AF); }

    private Set<Breakend> findRescuedBreakends(final LinkStore linkStore, final FilterCache filterCache, boolean rescueShortSVs)
    {
        mLinkRescue.getRescueInfo().clear();
        mLinkRescue.findRescuedBreakends(linkStore, filterCache, rescueShortSVs);
        return mLinkRescue.getRescueInfo().keySet();
    }

    private Set<Breakend> findRescuedDsbLineInsertions(final LinkStore linkStore, final FilterCache filterCache, double minQual)
    {
        mLinkRescue.getRescueInfo().clear();
        mLinkRescue.findRescuedDsbLineInsertions(linkStore, filterCache, minQual);
        return mLinkRescue.getRescueInfo().keySet();
    }

    @Test
    public void testVariantRescue()
    {
        SvData var1 = mGripss.createDel(CHR_1, 100, 200, null, null);
        SvData var2 = mGripss.createDel(CHR_1, 300, 400, null, null);
        SvData var3 = mGripss.createDel(CHR_1, 500, 600, null, null);
        SvData var4 = mGripss.createDel(CHR_1, 700, 800, null, null);

        setFiltered(var1.breakendStart());
        setFiltered(var1.breakendEnd());
        setFiltered(var2.breakendStart());
        setFiltered(var2.breakendEnd());
        setFiltered(var3.breakendStart());
        setFiltered(var3.breakendEnd());

        LinkStore linkStore = new LinkStore();

        linkStore.addLinks("1", var1.breakendStart(), var2.breakendStart());
        linkStore.addLinks("2", var2.breakendEnd(), var3.breakendEnd());
        linkStore.addLinks("3", var3.breakendStart(), var4.breakendStart());

        assertTrue(findRescuedBreakends(linkStore, mGripss.FilterCache, false).isEmpty());

        Set<Breakend> rescues = findRescuedBreakends(linkStore, mGripss.FilterCache, true);
        assertTrue(rescues.contains(var1.breakendStart()));
        assertTrue(rescues.contains(var1.breakendEnd()));
        assertTrue(rescues.contains(var2.breakendStart()));
        assertTrue(rescues.contains(var2.breakendEnd()));
        assertTrue(rescues.contains(var3.breakendStart()));
        assertTrue(rescues.contains(var3.breakendEnd()));

        // test can handle a loop of links
        linkStore.clear();
        linkStore.addLinks("1", var1.breakendEnd(), var2.breakendStart());
        linkStore.addLinks("2", var2.breakendEnd(), var3.breakendStart());
        linkStore.addLinks("3", var3.breakendEnd(), var4.breakendStart());
        linkStore.addLinks("4", var4.breakendEnd(), var1.breakendStart());

        rescues = findRescuedBreakends(linkStore, mGripss.FilterCache, true);
        assertTrue(rescues.contains(var1.breakendStart()));
        assertTrue(rescues.contains(var1.breakendEnd()));
        assertTrue(rescues.contains(var2.breakendStart()));
        assertTrue(rescues.contains(var2.breakendEnd()));
        assertTrue(rescues.contains(var3.breakendStart()));
        assertTrue(rescues.contains(var3.breakendEnd()));
    }

    @Test
    public void testRescueLineInsertions()
    {
        // rescue low qual SVs in a deletion bridge
        Map<String,Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(QUAL, 200);

        SvData var1 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 100, 20000, POS_ORIENT, NEG_ORIENT, LINE_INSERT_SEQ_A,
                mGripss.GenotypeIds, null, null, tumorAttributes);

        SvData var2 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 2010, 50000, POS_ORIENT, NEG_ORIENT, LINE_INSERT_SEQ_T,
                mGripss.GenotypeIds, null, null, tumorAttributes);

        assertTrue(var1.breakendEnd().IsLineInsertion);
        assertTrue(var2.breakendStart().IsLineInsertion);

        setFiltered(var1.breakendStart());
        setFiltered(var1.breakendEnd());
        setFiltered(var2.breakendStart());
        setFiltered(var2.breakendEnd());

        LinkStore linkStore = new LinkStore();

        linkStore.addLinks("dsb1", var1.breakendEnd(), var2.breakendStart());

        Set<Breakend> rescues = findRescuedDsbLineInsertions(linkStore, mGripss.FilterCache, 300);
        assertTrue(rescues.contains(var1.breakendEnd()));
        assertTrue(rescues.contains(var2.breakendStart()));

    }
}
