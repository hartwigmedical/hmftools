package com.hartwig.hmftools.sage.select;

import static com.hartwig.hmftools.sage.common.TestUtils.region;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

public class PanelSelectorTest
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
        assertEquals(ReadPanelStatus.WITHIN_PANEL, mPanelSelector.panelStatus(995));
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
}
