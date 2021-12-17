package com.hartwig.hmftools.common.genome.region;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class GenomeRegionsTest {

    private static final String CHROM = "1";

    @Test
    public void testAddPosition() {
        final GenomeRegions victim = new GenomeRegions(CHROM, 1);
        final List<Integer> positions = Lists.newArrayList(1, 1, 2, 3, 5, 6, 7, 9, 10, 10);
        Collections.shuffle(positions);

        for (Integer position : positions) {
            victim.addPosition(position);
        }

        final List<GenomeRegion> regions = victim.build();
        Assert.assertEquals(3, regions.size());
        assertRegion(regions.get(0), 1, 3);
        assertRegion(regions.get(1), 5, 7);
        assertRegion(regions.get(2), 9, 10);
    }

    @Test
    public void testInsideMinGapInFront() {
        final GenomeRegions victim = new GenomeRegions(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(4900);
        victim.addPosition(4950);

        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 4900, 5000);
    }

    @Test
    public void testOutsideMinGapInFront() {
        final GenomeRegions victim = new GenomeRegions(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(4899);

        Assert.assertEquals(2, victim.build().size());
        assertRegion(victim.build().get(0), 4899, 4899);
        assertRegion(victim.build().get(1), 5000, 5000);
    }

    @Test
    public void testJoinMinGapInFront() {
        final GenomeRegions victim = new GenomeRegions(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(4800);
        victim.addPosition(4900);

        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 4800, 5000);
    }

    @Test
    public void testInsideMinGapBehind() {
        final GenomeRegions victim = new GenomeRegions(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(5100);
        victim.addPosition(5050);

        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 5100);
    }

    @Test
    public void testOutsideMinGapBehind() {
        final GenomeRegions victim = new GenomeRegions(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(5101);

        Assert.assertEquals(2, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 5000);
        assertRegion(victim.build().get(1), 5101, 5101);
    }

    @Test
    public void testJoinMinGapBehind() {
        final GenomeRegions victim = new GenomeRegions(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(5200);
        victim.addPosition(5100);

        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 5200);
    }

    @Test
    public void testJoinRegionsWithMinGap() {
        GenomeRegions victim = new GenomeRegions(CHROM, 100);

        victim.addRegion(5000, 6000);
        victim.addRegion(6100, 7000);
        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 7000);

        victim = new GenomeRegions(CHROM, 100);

        victim.addRegion(5000, 6000);
        victim.addRegion(6101, 7000);
        Assert.assertEquals(2, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 6000);
        assertRegion(victim.build().get(1), 6101, 7000);
    }

    @Test
    public void testJoinRegions() {
        GenomeRegions victim = new GenomeRegions(CHROM, 1);

        victim.addRegion(5000, 6000);
        victim.addRegion(5500, 6500);
        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 6500);

        victim.addRegion(6501, 7000);
        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 7000);

        victim.addRegion(7002, 8000);
        Assert.assertEquals(2, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 7000);
        assertRegion(victim.build().get(1), 7002, 8000);
    }

    @Test
    public void testJoinRegionsBoundary() {
        GenomeRegions victim = new GenomeRegions(CHROM);

        victim.addRegion(5000, 5500);
        victim.addRegion(6000, 6500);
        Assert.assertEquals(2, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 5500);
        assertRegion(victim.build().get(1), 6000, 6500);

        victim.addRegion(5501, 5999);
        Assert.assertEquals(1, victim.build().size());
        assertRegion(victim.build().get(0), 5000, 6500);
    }

    private void assertRegion(@NotNull final GenomeRegion region, int expectedStart, int expectedEnd) {
        Assert.assertEquals(expectedStart, region.start());
        Assert.assertEquals(expectedEnd, region.end());
    }
}
