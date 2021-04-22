package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GCAccumulatorTest {

    private static final int WINDOW_SIZE = 1000;
    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1e-10;

    @Test
    public void testStandard() {
        final GenomeRegion region = GenomeRegions.create(CHROMOSOME, 1001, 10000);
        final GCAccumulator victim = new GCAccumulator(region);
        final GenomeRegionSelector<GCProfile> selector = selector(profile(1001, 0.30), profile(4001, 0.31), profile(9001, 0.35));
        selector.select(region, victim);
        assertEquals(0.32, victim.averageGCContent(), EPSILON);
    }

    @Test
    public void testExcludeUnmappable() {
        final GenomeRegion region = GenomeRegions.create(CHROMOSOME, 1001, 3000);
        final GCAccumulator victim = new GCAccumulator(region);
        final GenomeRegionSelector<GCProfile> selector = selector(profile(1001, 0.90), unmappableProfile(2001, 0.91));
        selector.select(region, victim);
        assertEquals(0.90, victim.averageGCContent(), EPSILON);
    }

    @Test
    public void testExcludeOverlapAtStart() {
        final GenomeRegion region = GenomeRegions.create(CHROMOSOME, 1100, 3000);
        final GCAccumulator victim = new GCAccumulator(region);
        final GenomeRegionSelector<GCProfile> selector = selector(profile(1001, 0.90), profile(2001, 0.91));
        selector.select(region, victim);
        assertEquals(0.91, victim.averageGCContent(), EPSILON);
    }

    @Test
    public void testExcludeOverlapFromEnd() {
        final GenomeRegion region = GenomeRegions.create(CHROMOSOME, 1001, 2999);
        final GCAccumulator victim = new GCAccumulator(region);
        final GenomeRegionSelector<GCProfile> selector = selector(profile(1001, 0.90), profile(2001, 0.91));
        selector.select(region, victim);
        assertEquals(0.90, victim.averageGCContent(), EPSILON);
    }

    @Test
    public void testExcludeFromBothEnds() {
        final GenomeRegion region = GenomeRegions.create(CHROMOSOME, 1100, 2999);
        final GCAccumulator victim = new GCAccumulator(region);
        final GenomeRegionSelector<GCProfile> selector = selector(profile(1001, 0.90), profile(2001, 0.91));
        selector.select(region, victim);
        assertEquals(0, victim.averageGCContent(), EPSILON);
    }

    @NotNull
    private static GenomeRegionSelector<GCProfile> selector(GCProfile... profiles) {
        List<GCProfile> list = Lists.newArrayList(profiles);
        return GenomeRegionSelectorFactory.create(list);
    }

    @NotNull
    private static GCProfile unmappableProfile(long start, double gcContent) {
        return create(start, GCProfile.MIN_MAPPABLE_PERCENTAGE - 0.1, gcContent);
    }

    @NotNull
    private static GCProfile profile(long start, double gcContent) {
        return create(start, GCProfile.MIN_MAPPABLE_PERCENTAGE, gcContent);
    }

    @NotNull
    private static GCProfile create(long start, double mappability, double gcContent) {
        return ImmutableGCProfile.builder()
                .chromosome(CHROMOSOME)
                .start(start)
                .end(start + WINDOW_SIZE - 1)
                .mappablePercentage(mappability)
                .gcContent(gcContent)
                .nonNPercentage(1)
                .build();
    }
}
