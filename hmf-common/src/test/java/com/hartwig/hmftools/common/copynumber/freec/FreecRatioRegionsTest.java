package com.hartwig.hmftools.common.copynumber.freec;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FreecRatioRegionsTest {

    @Test
    public void testSingleRegion() {
        final FreecRatio ratio = create("1", 1001, 0.5, 1);
        final List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(Lists.newArrayList(ratio));
        assertEquals(1, regions.size());
        assertRegion(regions.get(0), 1001, 2000);
    }

    @Test
    public void testMultipleSimilarRegions() {
        final FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        final FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        final FreecRatio ratio3 = create("1", 3001, 0.5, 1);
        final List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(
                Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(1, regions.size());
        assertRegion(regions.get(0), 1001, 4000);
    }

    @Test
    public void testMultipleRegionsDifferentBAF() {
        final FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        final FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        final FreecRatio ratio3 = create("1", 3001, 0.7, 1);
        final List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(
                Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(2, regions.size());
        assertRegion(regions.get(0), 1001, 3000);
        assertRegion(regions.get(1), 3001, 4000);
    }

    @Test
    public void testMultipleRegionsDifferentMedianRatio() {
        FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        FreecRatio ratio3 = create("1", 3001, 0.5, 1.1);
        List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(
                Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(2, regions.size());
        assertRegion(regions.get(0), 1001, 3000);
        assertRegion(regions.get(1), 3001, 4000);
    }

    @Test
    public void testMultipleRegionsDifferentChromosome() {
        FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        FreecRatio ratio3 = create("2", 1001, 0.5, 1);
        List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(
                Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(2, regions.size());
        assertRegion(regions.get(0), 1001, 3000);
        assertRegion(regions.get(1), 1001, 2000);
    }

    @Test
    public void testMultipleRegionsMissingSections() {
        final FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        final FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        final FreecRatio ratio3 = create("1", 4001, 0.5, 1);
        final List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(
                Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(1, regions.size());
        assertRegion(regions.get(0), 1001, 5000);
    }

    @Test
    public void testMultipleRegionsWithBreaks() {
        final FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        final FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        final FreecRatio ratio3 = create("1", 4001, 0.7, 1);
        final List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(
                Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(2, regions.size());
        assertRegion(regions.get(0), 1001, 3000);
        assertRegion(regions.get(1), 4001, 5000);
    }

    private static void assertRegion(@NotNull final GenomeRegion region, final long expectedStart,
            final long expectedEnd) {
        assertEquals(expectedStart, region.start());
        assertEquals(expectedEnd, region.end());
    }

    @NotNull
    private static FreecRatio create(@NotNull final String chromosome, final long position, final double estimatedBaf,
            final double medianRatio) {
        return ImmutableFreecRatio.builder().chromosome(chromosome).position(position).estimatedBAF(
                estimatedBaf).ratio(medianRatio).medianRatio(medianRatio).build();
    }
}
