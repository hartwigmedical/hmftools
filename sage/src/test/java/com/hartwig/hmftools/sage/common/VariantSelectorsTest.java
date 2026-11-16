package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.common.TestUtils.region;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.select.HotspotCache;
import com.hartwig.hmftools.sage.select.PanelSelector;
import com.hartwig.hmftools.sage.select.ReadPanelStatus;

import org.junit.Assert;
import org.junit.Test;

public class VariantSelectorsTest
{
    private final List<BaseRegion> panel =
            Lists.newArrayList(region(995, 995), region(998, 1102), region(1995, 1995), region(1998, 2102));

    private final PanelSelector mPanelSelector = new PanelSelector(panel);

    @Test
    public void testOverlap()
    {
        assertTrue(mPanelSelector.inPanel(995, 995));
        assertTrue(mPanelSelector.inPanel(994, 996));
        assertTrue(mPanelSelector.inPanel(1102, 1102));
        assertTrue(mPanelSelector.inPanel(1000, 1002));
        assertTrue(mPanelSelector.inPanel(996, 1000));
        assertFalse(mPanelSelector.inPanel(1, 994));
        assertTrue(mPanelSelector.inPanel(998, 998));
        assertFalse(mPanelSelector.inPanel(996, 997));
        assertTrue(mPanelSelector.inPanel(994, 1000));
        assertFalse(mPanelSelector.inPanel(2200, 10000));
    }

    @Test
    public void testOutOfOrder()
    {
        testOverlap();
        testOverlap();
    }

    @Test
    public void testPanelStatus()
    {
        Assert.assertEquals(ReadPanelStatus.WITHIN_PANEL, mPanelSelector.panelStatus(995));
        assertEquals(ReadPanelStatus.OUTSIDE_PANEL, mPanelSelector.panelStatus(996));
        assertEquals(ReadPanelStatus.OUTSIDE_PANEL, mPanelSelector.panelStatus(994));
        assertEquals(ReadPanelStatus.OUTSIDE_PANEL, mPanelSelector.panelStatus(2103));
        assertEquals(ReadPanelStatus.WITHIN_PANEL, mPanelSelector.panelStatus(2100));

        assertEquals(ReadPanelStatus.MIXED, mPanelSelector.panelStatus(995, 996));
        assertEquals(ReadPanelStatus.MIXED, mPanelSelector.panelStatus(994, 995));
        assertEquals(ReadPanelStatus.MIXED, mPanelSelector.panelStatus(994, 996)); // straddling
        assertEquals(ReadPanelStatus.MIXED, mPanelSelector.panelStatus(1995, 2000));
        assertEquals(ReadPanelStatus.MIXED, mPanelSelector.panelStatus(1500, 2500)); // straddling
    }

    @Test
    public void testHotspotSelection()
    {
        String ref = "A";
        List<SimpleVariant> hotspots = List.of(
                new SimpleVariant(CHR_1, 100, ref, "C"),
                new SimpleVariant(CHR_1, 101, ref, "C"),
                new SimpleVariant(CHR_1, 101, ref, "G"),
                new SimpleVariant(CHR_1, 101, ref, "T"),
                new SimpleVariant(CHR_1, 102, ref, "C"),
                new SimpleVariant(CHR_1, 103, ref, "C"),
                new SimpleVariant(CHR_1, 103, ref, "G"),
                new SimpleVariant(CHR_1, 103, ref, "T"),
                new SimpleVariant(CHR_1, 104, ref, "C"));

        HotspotCache hotspotCache = new HotspotCache(hotspots);

        SimpleVariant testVariant = new SimpleVariant(CHR_1, 99, ref, "C");
        assertFalse(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 105, ref, "C");
        assertFalse(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 100, ref, "G");
        assertFalse(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 100, ref, "C");
        assertTrue(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 101, ref, "C");
        assertTrue(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 101, ref, "T");
        assertTrue(hotspotCache.isHotspot(testVariant));

        // moving backwards through matches
        testVariant = new SimpleVariant(CHR_1, 101, ref, "G");
        assertTrue(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 104, ref, "G");
        assertFalse(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 104, ref, "C");
        assertTrue(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 103, ref, "C");
        assertTrue(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 103, ref, "T");
        assertTrue(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 103, ref, "G");
        assertTrue(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 102, ref, "G");
        assertFalse(hotspotCache.isHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 102, ref, "C");
        assertTrue(hotspotCache.isHotspot(testVariant));
    }

    @Test
    public void testHotspotProximity()
    {
        String ref = "A";
        String alt = "C";

        List<SimpleVariant> hotspots = List.of(
                new SimpleVariant(CHR_1, 100, ref, alt),
                new SimpleVariant(CHR_1, 110, ref, alt),
                new SimpleVariant(CHR_1, 110, ref, alt),
                new SimpleVariant(CHR_1, 122, ref, alt),
                new SimpleVariant(CHR_1, 123, ref, alt),
                new SimpleVariant(CHR_1, 123, ref, alt),
                new SimpleVariant(CHR_1, 123, "AAA", alt),
                new SimpleVariant(CHR_1, 140, ref, alt));

        HotspotCache hotspotCache = new HotspotCache(hotspots);

        SimpleVariant testVariant = new SimpleVariant(CHR_1, 94, ref, alt);
        assertFalse(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 95, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 105, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 116, ref, alt);
        assertFalse(hotspotCache.isNearHotspot(testVariant));

        // test navigation
        testVariant = new SimpleVariant(CHR_1, 124, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 120, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 135, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 145, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 146, ref, alt);
        assertFalse(hotspotCache.isNearHotspot(testVariant));

        hotspotCache.setIndex(0);

        testVariant = new SimpleVariant(CHR_1, 124, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        hotspotCache.setIndex(hotspots.size() - 1);

        testVariant = new SimpleVariant(CHR_1, 119, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        // check indels
        testVariant = new SimpleVariant(CHR_1, 130, ref, alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));

        testVariant = new SimpleVariant(CHR_1, 133, "AAA", alt);
        assertTrue(hotspotCache.isNearHotspot(testVariant));
    }
}
