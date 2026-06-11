package com.hartwig.hmftools.redux.unmap;

import static com.hartwig.hmftools.redux.ReduxConstants.UNMAP_MIN_HIGH_DEPTH;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.mappability.UnmappingRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class ReadUnmapperMergeTest
{
    private static final String CHR_1 = "1";

    private static Map<String,List<ChrBaseRegion>> rnaRegions(final ChrBaseRegion... regions)
    {
        final Map<String,List<ChrBaseRegion>> map = Maps.newHashMap();
        for(ChrBaseRegion region : regions)
            map.computeIfAbsent(region.Chromosome, x -> Lists.newArrayList()).add(region);
        return map;
    }

    @Test
    public void addsRegionToEmptyChromosome()
    {
        final ReadUnmapper unmapper = new ReadUnmapper(Maps.newHashMap());
        unmapper.mergeAlwaysUnmapRegions(rnaRegions(new ChrBaseRegion(CHR_1, 500, 700)));

        final List<UnmappingRegion> regions = unmapper.getRegions(CHR_1);
        assertEquals(1, regions.size());
        assertEquals(500, regions.get(0).start());
        assertEquals(700, regions.get(0).end());
        assertEquals(UNMAP_MIN_HIGH_DEPTH, regions.get(0).maxDepth());
    }

    @Test
    public void overlappingMappabilityRegionIsAbsorbedAsHighDepth()
    {
        final ReadUnmapper unmapper = new ReadUnmapper(Maps.newHashMap());
        // existing low-depth mappability region inside the curated zone
        unmapper.addRegion(CHR_1, new UnmappingRegion(550, 600, 10));

        unmapper.mergeAlwaysUnmapRegions(rnaRegions(new ChrBaseRegion(CHR_1, 500, 700)));

        final List<UnmappingRegion> regions = unmapper.getRegions(CHR_1);
        assertEquals(1, regions.size());
        assertEquals(500, regions.get(0).start());
        assertEquals(700, regions.get(0).end());
        // high-depth wins so the merged region unmaps unconditionally
        assertEquals(UNMAP_MIN_HIGH_DEPTH, regions.get(0).maxDepth());
    }

    @Test
    public void nonOverlappingRegionsKeptSortedAndDistinct()
    {
        final ReadUnmapper unmapper = new ReadUnmapper(Maps.newHashMap());
        unmapper.addRegion(CHR_1, new UnmappingRegion(100, 200, 10));
        unmapper.addRegion(CHR_1, new UnmappingRegion(900, 1000, 20));

        // curated zone between the two existing regions
        unmapper.mergeAlwaysUnmapRegions(rnaRegions(new ChrBaseRegion(CHR_1, 500, 600)));

        final List<UnmappingRegion> regions = unmapper.getRegions(CHR_1);
        assertEquals(3, regions.size());
        assertEquals(100, regions.get(0).start());
        assertEquals(500, regions.get(1).start());
        assertEquals(UNMAP_MIN_HIGH_DEPTH, regions.get(1).maxDepth());
        assertEquals(900, regions.get(2).start());
    }
}
