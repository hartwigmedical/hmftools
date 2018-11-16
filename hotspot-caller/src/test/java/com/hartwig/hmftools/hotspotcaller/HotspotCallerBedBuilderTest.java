package com.hartwig.hmftools.hotspotcaller;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class HotspotCallerBedBuilderTest {

    private static final String CHROM = "1";

    @Test
    public void testAddPosition() {
        final List<Long> positions = Lists.newArrayList(1L, 1L, 2L, 3L, 5L, 6L, 7L, 9L, 10L, 10L);
        Collections.shuffle(positions);

        List<GenomeRegion> regions = Lists.newArrayList();
        for (Long position : positions) {
            regions = HotspotCallerBedBuilder.addPosition(CHROM, position, regions);
        }

        Assert.assertEquals(3, regions.size());
        assertRegion(regions.get(0), 1, 3);
        assertRegion(regions.get(1), 5, 7);
        assertRegion(regions.get(2), 9, 10);
    }

    @Test
    public void testAddMNV() {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome(CHROM).alt("GAT").ref("TAC").position(10).build();
        final List<GenomeRegion> regions = HotspotCallerBedBuilder.addVariantHotspot(hotspot, Lists.newArrayList());
        Assert.assertEquals(1, regions.size());
        assertRegion(regions.get(0), 10, 12);
    }

    @Test
    public void testAddIndel() {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome(CHROM).alt("G").ref("TAC").position(10).build();
        final List<GenomeRegion> regions = HotspotCallerBedBuilder.addVariantHotspot(hotspot, Lists.newArrayList());
        Assert.assertEquals(1, regions.size());
        assertRegion(regions.get(0), 10, 10);
    }

    private void assertRegion(@NotNull final GenomeRegion region, long expectedStart, long expectedEnd) {
        Assert.assertEquals(expectedStart, region.start());
        Assert.assertEquals(expectedEnd, region.end());
    }

}
