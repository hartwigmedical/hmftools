package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class TargetDepthTest
{
    @Test
    public void testTargetDepthTracker()
    {
        int readLength = 100;
        int targetDepth = 5;
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1000, 2000);

        TargetDepthTracker targetDepthTracker = new TargetDepthTracker(region, true, targetDepth, 100);

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
        targetDepthTracker = new TargetDepthTracker(region, true, targetDepth, 100);

        for(int i = 0; i < targetDepth; ++i)
        {
            assertTrue(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
        }

        assertFalse(targetDepthTracker.processReadDepth(false, position, position + readLength - 1));
    }
}
