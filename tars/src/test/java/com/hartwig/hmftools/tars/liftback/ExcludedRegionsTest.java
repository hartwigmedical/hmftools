package com.hartwig.hmftools.tars.liftback;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

// excludes() overlap contract
public class ExcludedRegionsTest
{
    private static final String CHR = "1";

    private static ExcludedRegions regions(final int start, final int end)
    {
        Map<String, List<BaseRegion>> map = new HashMap<>();
        map.put(CHR, new ArrayList<>(List.of(new BaseRegion(start, end))));
        return new ExcludedRegions(map);
    }

    @Test
    public void testExcludes()
    {
        // span inside the region
        assertTrue(regions(1000, 2000).excludes(CHR, 1500, 1600));

        // span outside the region
        assertFalse(regions(1000, 2000).excludes(CHR, 2500, 2600));

        // different chromosome
        assertFalse(regions(1000, 2000).excludes("2", 1500, 1600));

        // query end == region start: inclusive overlap fires
        assertTrue(regions(2000, 3000).excludes(CHR, 1900, 2000));

        // query start == region end: inclusive overlap fires
        assertTrue(regions(1000, 2000).excludes(CHR, 2000, 2100));
    }
}
