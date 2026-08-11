package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class TargetDepthTest
{
    @Test
    public void testTargetDepthTrackerWindowCoverage()
    {
        int readLength = 100;
        int targetDepth = 5;
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1000, 2000);

        TargetDepthTracker targetDepthTracker = new TargetDepthTracker(
                region, true, targetDepth, 100, true);

        for(int i = 0; i < 10; ++i)
        {
            assertTrue(targetDepthTracker.processReadDepth(false, 900, 950));
        }

        int position = 950;
        for(int i = 0; i < targetDepth; ++i)
        {
            assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
            position += 10;
        }

        assertFalse(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));

        position = 1051;
        assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));

        assertFalse(targetDepthTracker.processReadDepth(false, position, position + 20));

        position = 1150;

        for(int i = 0; i < 4; ++i)
        {
            assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
        }

        // try again with a single-base region
        region = new ChrBaseRegion(CHR_1, 1000, 1000);
        position = 980;
        targetDepthTracker = new TargetDepthTracker(
                region, true, targetDepth, 100, true);

        for(int i = 0; i < targetDepth; ++i)
        {
            assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
        }

        assertFalse(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
    }

    @Test
    public void testTargetDepthTrackerRollingWindow()
    {
        int readLength = 100;
        int targetDepth = 20;
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1000, 2000);

        TargetDepthTracker targetDepthTracker = new TargetDepthTracker(
                region, true, targetDepth, 100, false);

        assertEquals(50, targetDepthTracker.windowSize());
        assertEquals(10, targetDepthTracker.baseMaxReadsPerWindow());

        // first reads within the same window up to the limit
        int position = 1000;
        int readCount = targetDepth / 2;
        for(int i = 0; i < readCount; ++i)
        {
            assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
            position += 4;
        }

        assertEquals(readCount, targetDepthTracker.windowReadCount());

        assertFalse(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));

        position = 1051;
        assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));

        assertEquals(1, targetDepthTracker.windowReadCount());

        // move onto a new window, carrying over spare capacity
        position = 1110;
        assertTrue(targetDepthTracker.processReadDepth(false, position, position + 20));

        assertEquals(13, targetDepthTracker.maxReadsPerWindow());


        // try again with a small region
        region = new ChrBaseRegion(CHR_1, 1000, 1030);
        targetDepthTracker = new TargetDepthTracker(
                region, true, targetDepth, 100, false);

        position = 980;
        for(int i = 0; i < targetDepth; ++i)
        {
            assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
            position += 2;
        }

        assertFalse(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
    }
}
