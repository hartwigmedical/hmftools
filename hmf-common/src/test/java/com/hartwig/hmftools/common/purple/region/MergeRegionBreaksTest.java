package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class MergeRegionBreaksTest {

    private static final String CHROM1 = "1";
    private static final double EPSILON = 1e-10;

    @Test
    public void noMerge() {
        ConsolidatedRegion region1 = create(CHROM1, 1, 10000, 0.3, 4);
        ConsolidatedRegion region2 = create(CHROM1, 10001, 11001, 0.31, 4.3);
        ConsolidatedRegion region3 = create(CHROM1, 11002, 20000, 0.32, 5);

        List<ConsolidatedRegion> regions = MergeRegionBreaks.merge(Lists.newArrayList(region1, region2, region3));
        assertEquals(3, regions.size());
        assertRegion(region1, regions.get(0));
        assertRegion(region2, regions.get(1));
        assertRegion(region3, regions.get(2));
    }

    @Test
    public void mergeDownLeft() {
        ConsolidatedRegion region0 = create(CHROM1, 1, 5000, 0.3, 4);
        ConsolidatedRegion region1 = create(CHROM1, 5001, 10000, 0.3, 4);
        ConsolidatedRegion region2 = create(CHROM1, 10001, 11000, 0.31, 4.3);
        ConsolidatedRegion region3 = create(CHROM1, 11001, 20000, 0.32, 5);

        List<ConsolidatedRegion> regions = MergeRegionBreaks.merge(Lists.newArrayList(region0, region1, region2, region3));
        assertEquals(3, regions.size());
        assertRegion(region0, regions.get(0));
        assertRegion(create(CHROM1, 5001, 11000, 0.3, 4), regions.get(1));
        assertRegion(region3, regions.get(2));
    }

    @Test
    public void mergeDownRight() {
        ConsolidatedRegion region1 = create(CHROM1, 1, 10000, 0.3, 5);
        ConsolidatedRegion region2 = create(CHROM1, 10001, 11000, 0.31, 4.3);
        ConsolidatedRegion region3 = create(CHROM1, 11001, 20000, 0.32, 4);
        ConsolidatedRegion region4 = create(CHROM1, 20001, 30000, 0.32, 4.1);

        List<ConsolidatedRegion> regions = MergeRegionBreaks.merge(Lists.newArrayList(region1, region2, region3, region4));
        assertEquals(3, regions.size());
        assertRegion(region1, regions.get(0));
        assertRegion(create(CHROM1, 10001, 20000, 0.32, 4), regions.get(1));
        assertRegion(region4, regions.get(2));
    }

    @Test
    public void mergeUpLeft() {
        ConsolidatedRegion region1 = create(CHROM1, 1, 10000, 0.3, 5);
        ConsolidatedRegion region2 = create(CHROM1, 10001, 11000, 0.31, 4.7);
        ConsolidatedRegion region3 = create(CHROM1, 11001, 20000, 0.32, 4);

        List<ConsolidatedRegion> regions = MergeRegionBreaks.merge(Lists.newArrayList(region1, region2, region3));
        assertEquals(2, regions.size());
        assertRegion(create(CHROM1, 1, 11000, 0.3, 5), regions.get(0));
        assertRegion(region3, regions.get(1));
    }

    @Test
    public void mergeUpRight() {
        ConsolidatedRegion region1 = create(CHROM1, 1, 10000, 0.3, 4);
        ConsolidatedRegion region2 = create(CHROM1, 10001, 11000, 0.31, 4.7);
        ConsolidatedRegion region3 = create(CHROM1, 11001, 20000, 0.32, 5);

        List<ConsolidatedRegion> regions = MergeRegionBreaks.merge(Lists.newArrayList(region1, region2, region3));
        assertEquals(2, regions.size());
        assertRegion(region1, regions.get(0));
        assertRegion(create(CHROM1, 10001, 20000, 0.32, 5), regions.get(1));
    }

    private void assertRegion(ConsolidatedRegion expected, ConsolidatedRegion victim) {
        assertEquals(expected.chromosome(), victim.chromosome());
        assertEquals(expected.start(), victim.start());
        assertEquals(expected.end(), victim.end());
        assertEquals(expected.averageObservedBAF(), victim.averageObservedBAF(), EPSILON);
        assertEquals(expected.averageTumorCopyNumber(), victim.averageTumorCopyNumber(), EPSILON);
    }


    private ConsolidatedRegion create(String chromosome, long start, long end, double baf, double ratio) {
        return ImmutableConsolidatedRegion.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageObservedBAF(baf)
                .averageTumorCopyNumber(ratio)
                .build();
    }

}
