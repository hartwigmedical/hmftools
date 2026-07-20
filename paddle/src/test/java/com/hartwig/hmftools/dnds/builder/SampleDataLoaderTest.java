package com.hartwig.hmftools.dnds.builder;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

public class SampleDataLoaderTest
{
    private static final String CHR_1 = "1";

    private static Map<String, List<ChrBaseRegion>> buildTargetRegions()
    {
        List<ChrBaseRegion> chr1Regions = Lists.newArrayList(
                new ChrBaseRegion(CHR_1, 100, 200),
                new ChrBaseRegion(CHR_1, 500, 600),
                new ChrBaseRegion(CHR_1, 1000, 1100),
                new ChrBaseRegion(CHR_1, 2000, 2100),
                new ChrBaseRegion(CHR_1, 3000, 3100));

        Map<String, List<ChrBaseRegion>> regionsByChromosome = Maps.newHashMap();
        regionsByChromosome.put(CHR_1, chr1Regions);
        return regionsByChromosome;
    }

    @Test
    public void testInTargetRegions()
    {
        Map<String, List<ChrBaseRegion>> targetRegions = buildTargetRegions();

        // hit, no narrowing required (middle region on first mid-point check)
        assertTrue(SampleDataLoader.inTargetRegions(targetRegions, CHR_1, 1050));

        // hit, requires narrowing to the lower half
        assertTrue(SampleDataLoader.inTargetRegions(targetRegions, CHR_1, 550));

        // miss, gap between regions
        assertFalse(SampleDataLoader.inTargetRegions(targetRegions, CHR_1, 700));

        // hits on region boundaries
        assertTrue(SampleDataLoader.inTargetRegions(targetRegions, CHR_1, 100));
        assertTrue(SampleDataLoader.inTargetRegions(targetRegions, CHR_1, 200));

        // miss, before the first region
        assertFalse(SampleDataLoader.inTargetRegions(targetRegions, CHR_1, 50));

        // miss, after the last region
        assertFalse(SampleDataLoader.inTargetRegions(targetRegions, CHR_1, 3200));

        // miss, chromosome not present in map at all
        assertFalse(SampleDataLoader.inTargetRegions(targetRegions, "2", 1050));

        // hit, VCF contig carries a 'chr' prefix while target regions are keyed without one (GRCh38 vs BED convention)
        assertTrue(SampleDataLoader.inTargetRegions(targetRegions, "chr1", 1050));

        // miss, non-human contig (e.g. unplaced scaffold) must not throw
        assertFalse(SampleDataLoader.inTargetRegions(targetRegions, "chrUn_KI270742v1", 1050));
    }
}
