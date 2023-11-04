package com.hartwig.hmftools.common.genome.region;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

public class BaseRegionTest
{
    @Test
    public void testBinarySearch()
    {
        int startPos = 1000;
        int regionLength = 100;
        int regionGap = 1000;
        int regionCount = 20;

        List<BaseRegion> regions = Lists.newArrayList();

        for(int i = 0; i < regionCount; ++i)
        {
            int regionStart = startPos + i * regionGap;
            int regionEnd = regionStart + regionLength - 1;

            regions.add(new BaseRegion(regionStart, regionEnd));
        }

        // before the first region still returns 1
        assertEquals(0, BaseRegion.binarySearch(1, regions));

        // exact start region start match
        assertEquals(0, BaseRegion.binarySearch(1000, regions));

        // within start region
        assertEquals(0, BaseRegion.binarySearch(1050, regions));

        // exact region start match
        assertEquals(8, BaseRegion.binarySearch(9000, regions));

        // within a region in the middle
        assertEquals(9, BaseRegion.binarySearch(10050, regions));

        // in the middle in between
        assertEquals(9, BaseRegion.binarySearch(10500, regions));

        // last region
        assertEquals(19, BaseRegion.binarySearch(20050, regions));

        // beyond the last region
        assertEquals(19, BaseRegion.binarySearch(100000, regions));
    }
}
