package com.hartwig.hmftools.sage;

import java.util.List;

import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class SageHotspotBedBuilderTest {

    private static final String CHROM = "1";

    @Test
    public void testAddMNV() {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome(CHROM).alt("GAT").ref("TAC").position(10).build();
        final List<GenomeRegion> regions = SageHotspotBedBuilder.addVariantHotspot(hotspot);
        Assert.assertEquals(1, regions.size());
        assertRegion(regions.get(0), 10, 12);
    }

    @Test
    public void testAddIndel() {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome(CHROM).alt("G").ref("TAC").position(10).build();
        final List<GenomeRegion> regions = SageHotspotBedBuilder.addVariantHotspot(hotspot);
        Assert.assertEquals(1, regions.size());
        assertRegion(regions.get(0), 10, 10);
    }

    private void assertRegion(@NotNull final GenomeRegion region, long expectedStart, long expectedEnd) {
        Assert.assertEquals(expectedStart, region.start());
        Assert.assertEquals(expectedEnd, region.end());
    }

}
