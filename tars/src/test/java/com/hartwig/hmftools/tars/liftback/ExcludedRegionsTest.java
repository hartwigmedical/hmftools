package com.hartwig.hmftools.tars.liftback;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

// excludes() overlap contract, queried post-lift on lifted genomic coords.
public class ExcludedRegionsTest
{
    private static final String CHR = "1";

    private static ExcludedRegions regions(final int start, final int end)
    {
        Map<String, List<ChrBaseRegion>> map = new HashMap<>();
        map.put(CHR, new ArrayList<>(List.of(new ChrBaseRegion(CHR, start, end))));
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

        // end just touches the region start -- inclusive overlap must fire
        assertTrue(regions(2000, 3000).excludes(CHR, 1900, 2000));

        // start just touches the region end -- inclusive overlap must fire
        assertTrue(regions(1000, 2000).excludes(CHR, 2000, 2100));
    }
}
