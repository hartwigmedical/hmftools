package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.buildLinkAttributes;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_CIPOS;
import static com.hartwig.hmftools.gripss.links.DsbLinkFinder.findBreaks;

import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertNull;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.links.AssemblyLinks;
import com.hartwig.hmftools.gripss.links.LinkStore;

import org.junit.Test;

public class DsbLinksTest
{
    private final GripssTestApplication mGripss;

    public DsbLinksTest()
    {
        mGripss = new GripssTestApplication();
    }

    @Test
    public void testBasicDsbLinks()
    {
        // basic DSB link
        SvData var1 = mGripss.createInv(CHR_1, 1000, 2000, POS_ORIENT, null, null);
        SvData var2 = mGripss.createInv(CHR_1, 2010, 3000, NEG_ORIENT, null, null);

        LinkStore assemblyLinks = AssemblyLinks.buildAssembledLinks(Lists.newArrayList(var1, var2));
        assertTrue(assemblyLinks.getBreakendLinksMap().isEmpty());

        SvDataCache dataCache = new SvDataCache();

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2));

        Set<Breakend> duplicateBreakends = Sets.newHashSet();
        LinkStore dsbLinks = findBreaks(dataCache, assemblyLinks, duplicateBreakends);

        assertNull(dsbLinks.getBreakendLinks(var1.breakendStart()));
        assertNotNull(dsbLinks.getBreakendLinks(var1.breakendEnd()));
        assertNotNull(dsbLinks.getBreakendLinks(var2.breakendStart()));
        assertNull(dsbLinks.getBreakendLinks(var2.breakendEnd()));

        dsbLinks.clear();

        // breakends can overlap as long as not assembled into a TI
        var1 = mGripss.createInv(CHR_1, 1000, 2010, POS_ORIENT, null, null);
        var2 = mGripss.createInv(CHR_1, 2000, 3000, NEG_ORIENT, null, null);

        assemblyLinks = AssemblyLinks.buildAssembledLinks(Lists.newArrayList(var1, var2));
        assertTrue(assemblyLinks.getBreakendLinksMap().isEmpty());

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2));

        dsbLinks = findBreaks(dataCache, assemblyLinks, duplicateBreakends);

        assertNull(dsbLinks.getBreakendLinks(var1.breakendStart()));
        assertNotNull(dsbLinks.getBreakendLinks(var1.breakendEnd()));
        assertNotNull(dsbLinks.getBreakendLinks(var2.breakendStart()));
        assertNull(dsbLinks.getBreakendLinks(var2.breakendEnd()));

        // will not link if already assembled
        var1 = mGripss.createInv(CHR_1, 1000, 2010, POS_ORIENT, null, buildLinkAttributes("asm1", "325"));
        var2 = mGripss.createInv(CHR_1, 2000, 3000, NEG_ORIENT, buildLinkAttributes("asm1", "325"), null);

        assemblyLinks = AssemblyLinks.buildAssembledLinks(Lists.newArrayList(var1, var2));
        assertTrue(!assemblyLinks.getBreakendLinksMap().isEmpty());

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2));

        dsbLinks = findBreaks(dataCache, assemblyLinks, duplicateBreakends);
        assertTrue(dsbLinks.getBreakendLinksMap().isEmpty());

        // will not link if too many alternatives
        var1 = mGripss.createInv(CHR_1, 1000, 2010, POS_ORIENT, null, null);
        var2 = mGripss.createInv(CHR_1, 2000, 3000, NEG_ORIENT, null, null);
        SvData var3 = mGripss.createInv(CHR_1, 2030, 5000, NEG_ORIENT, null, null);
        SvData var4 = mGripss.createInv(CHR_1, 2055, 8000, POS_ORIENT, null, null);

        assemblyLinks = AssemblyLinks.buildAssembledLinks(Lists.newArrayList(var1, var2, var3, var4));
        assertTrue(assemblyLinks.getBreakendLinksMap().isEmpty());

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2, var3, var4));

        dsbLinks = findBreaks(dataCache, assemblyLinks, duplicateBreakends);
        assertTrue(dsbLinks.getBreakendLinksMap().isEmpty());
    }

    @Test
    public void testDsbMultipleBreakends()
    {
        // skips over breakend which is too far away but links to one which has a confidence interval that brings it close enough
        SvData var1 = mGripss.createInv(CHR_1, 110, 2000, POS_ORIENT, null, null);
        SvData var2 = mGripss.createInv(CHR_1, 220, 5000, NEG_ORIENT, null, null);

        Map<String,Object> attributes = Maps.newHashMap();
        attributes.put(VT_CIPOS, new int[] {-400, 0});
        SvData var3 = mGripss.createInv(CHR_1, 420, 10000, NEG_ORIENT, attributes, null);

        LinkStore assemblyLinks = AssemblyLinks.buildAssembledLinks(Lists.newArrayList(var1, var2, var3));
        assertTrue(assemblyLinks.getBreakendLinksMap().isEmpty());

        SvDataCache dataCache = new SvDataCache();

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2, var3));

        Set<Breakend> duplicateBreakends = Sets.newHashSet();
        LinkStore dsbLinks = findBreaks(dataCache, assemblyLinks, duplicateBreakends);

        assertNotNull(dsbLinks.getBreakendLinks(var1.breakendStart()));
        assertNotNull(dsbLinks.getBreakendLinks(var3.breakendStart()));
        assertNull(dsbLinks.getBreakendLinks(var2.breakendStart()));
    }

    @Test
    public void testDsbIgnoreDuplicates()
    {
        // skips over breakend which is too far away but links to one which has a confidence interval that brings it close enough
        SvData var1 = mGripss.createInv(CHR_1, 110, 2000, POS_ORIENT, null, null);
        SvData var2 = mGripss.createInv(CHR_1, 130, 5000, POS_ORIENT, null, null);
        SvData var3 = mGripss.createInv(CHR_1, 150, 8000, NEG_ORIENT, null, null);

        LinkStore assemblyLinks = AssemblyLinks.buildAssembledLinks(Lists.newArrayList(var1, var2, var3));
        assertTrue(assemblyLinks.getBreakendLinksMap().isEmpty());

        SvDataCache dataCache = new SvDataCache();

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2, var3));

        Set<Breakend> duplicateBreakends = Sets.newHashSet(var1.breakendStart());
        LinkStore dsbLinks = findBreaks(dataCache, assemblyLinks, duplicateBreakends);

        assertNull(dsbLinks.getBreakendLinks(var1.breakendStart()));
        assertNotNull(dsbLinks.getBreakendLinks(var2.breakendStart()));
        assertNotNull(dsbLinks.getBreakendLinks(var3.breakendStart()));
    }
}
