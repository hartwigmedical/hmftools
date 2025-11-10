package com.hartwig.hmftools.panelbuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class TargetedRangeTest
{
    @Test
    public void testFromStartAndLength()
    {
        assertEquals(new TargetedRange(30, 130), TargetedRange.fromStartAndLength(30, 100));
    }

    @Test
    public void testSingleBase()
    {
        assertEquals(new TargetedRange(32, 33), TargetedRange.singleBase(32));
    }

    @Test
    public void testWholeRegion()
    {
        assertEquals(new TargetedRange(0, 123), TargetedRange.wholeRegion(123));
    }

    @Test
    public void testFromRegions()
    {
        assertEquals(
                new TargetedRange(20, 101),
                TargetedRange.fromRegions(new ChrBaseRegion("1", 120, 230), new ChrBaseRegion("1", 100, 200)));
    }

    @Test
    public void testFromRegionsInvalid()
    {
        assertThrows(RuntimeException.class,
                () -> TargetedRange.fromRegions(new ChrBaseRegion("1", 100, 200), new ChrBaseRegion("2", 100, 200)));
        assertThrows(RuntimeException.class,
                () -> TargetedRange.fromRegions(new ChrBaseRegion("1", 100, 200), new ChrBaseRegion("1", 201, 300)));

    }
}
