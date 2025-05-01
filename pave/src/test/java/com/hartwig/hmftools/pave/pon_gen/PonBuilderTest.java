package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.junit.Test;

public class PonBuilderTest
{
    @Test
    public void testAddNewVariants()
    {
        VariantCache variantCache = new VariantCache();

        VariantPonData variant = variantCache.getOrCreateVariant(CHR_1, 100, "A", "C");
        assertEquals(1, variantCache.variantCount());

        variantCache.getOrCreateVariant(CHR_1, 100, "A", "C");
        assertEquals(1, variantCache.variantCount());

        variantCache.getOrCreateVariant(CHR_1, 100, "AA", "C");
        assertEquals(2, variantCache.variantCount());

        variantCache.getOrCreateVariant(CHR_1, 100, "A", "G");
        assertEquals(3, variantCache.variantCount());

        variantCache.getOrCreateVariant(CHR_1, 101, "A", "G");
        assertEquals(4, variantCache.variantCount());
        assertEquals(3, variantCache.lastVariantIndex());

        // will warn on variants out of order (ie earlier) but still handle this
        variantCache.getOrCreateVariant(CHR_1, 99, "A", "G");
        assertEquals(5, variantCache.variantCount());
        assertEquals(0, variantCache.lastVariantIndex());
    }

    @Test
    public void testHotspotCache()
    {
        VariantHotspot hotspot1 = ImmutableVariantHotspotImpl.builder().chromosome(CHR_1).position(100).ref("A").alt("C").build();
        VariantHotspot hotspot2 = ImmutableVariantHotspotImpl.builder().chromosome(CHR_1).position(100).ref("A").alt("AA").build();
        VariantHotspot hotspot3 = ImmutableVariantHotspotImpl.builder().chromosome(CHR_1).position(100).ref("A").alt("G").build();
        VariantHotspot hotspot4 = ImmutableVariantHotspotImpl.builder().chromosome(CHR_1).position(105).ref("A").alt("C").build();
        VariantHotspot hotspot5 = ImmutableVariantHotspotImpl.builder().chromosome(CHR_1).position(110).ref("A").alt("C").build();

        List<VariantHotspot> hotspots = Lists.newArrayList(hotspot1, hotspot2, hotspot3, hotspot4, hotspot5);

        HotspotRegionCache hotspotCache = new HotspotRegionCache(hotspots);

        assertFalse(hotspotCache.matchesHotspot(95, "A", "C"));

        assertTrue(hotspotCache.matchesHotspot(100, "A", "C"));
        assertFalse(hotspotCache.matchesHotspot(100, "A", "T"));

        assertTrue(hotspotCache.matchesHotspot(100, "A", "G"));
        assertTrue(hotspotCache.matchesHotspot(100, "A", "AA"));
        assertTrue(hotspotCache.matchesHotspot(100, "A", "C"));

        assertFalse(hotspotCache.matchesHotspot(102, "A", "C"));

        assertTrue(hotspotCache.matchesHotspot(105, "A", "C"));
        assertTrue(hotspotCache.matchesHotspot(110, "A", "C"));
        assertEquals(4, hotspotCache.currentIndex());

        // cannot go back
        assertFalse(hotspotCache.matchesHotspot(100, "A", "C"));
        assertFalse(hotspotCache.matchesHotspot(105, "A", "C"));

        hotspotCache.resetSearch();
        assertTrue(hotspotCache.matchesHotspot(100, "A", "C"));
        assertTrue(hotspotCache.matchesHotspot(100, "A", "AA"));
        assertEquals(0, hotspotCache.currentIndex());
    }

    @Test
    public void testExonicRegionCache()
    {
        List<BaseRegion> regions = Lists.newArrayList(
                new BaseRegion(100, 200),
                new BaseRegion(500, 1000),
                new BaseRegion(2000, 2200));

        TranscriptRegionCache transcriptRegionCache = new TranscriptRegionCache(regions);

        assertFalse(transcriptRegionCache.inOrNearExonicRegion(95));
        assertTrue(transcriptRegionCache.inOrNearExonicRegion(100));
        assertTrue(transcriptRegionCache.inOrNearExonicRegion(200));

        assertFalse(transcriptRegionCache.inOrNearExonicRegion(250));
        assertTrue(transcriptRegionCache.inOrNearExonicRegion(750));

        // cannot go back
        assertFalse(transcriptRegionCache.inOrNearExonicRegion(150));

        assertFalse(transcriptRegionCache.inOrNearExonicRegion(5000));
    }
}
