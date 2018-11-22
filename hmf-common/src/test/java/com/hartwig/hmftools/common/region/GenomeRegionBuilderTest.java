package com.hartwig.hmftools.common.region;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class GenomeRegionBuilderTest {

    private static final String CHROM = "1";

    @Test
    public void testAddPosition() {
        final GenomeRegionBuilder victim = new GenomeRegionBuilder(CHROM);
        final List<Long> positions = Lists.newArrayList(1L, 1L, 2L, 3L, 5L, 6L, 7L, 9L, 10L, 10L);
        Collections.shuffle(positions);

        for (Long position : positions) {
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

        final GenomeRegionBuilder victim = new GenomeRegionBuilder(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(4900);
        victim.addPosition(4950);

        Assert.assertEquals(1, victim.build().size());
        assertRegion( victim.build().get(0), 4900, 5000);
    }

    @Test
    public void testOutsideMinGapInFront() {

        final GenomeRegionBuilder victim = new GenomeRegionBuilder(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(4899);

        Assert.assertEquals(2, victim.build().size());
        assertRegion( victim.build().get(0), 4899, 4899);
        assertRegion( victim.build().get(1), 5000, 5000);
    }

    @Test
    public void testJoinMinGapInFront() {

        final GenomeRegionBuilder victim = new GenomeRegionBuilder(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(4800);
        victim.addPosition(4900);

        Assert.assertEquals(1, victim.build().size());
        assertRegion( victim.build().get(0), 4800, 5000);
    }

    @Test
    public void testInsideMinGapBehind() {

        final GenomeRegionBuilder victim = new GenomeRegionBuilder(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(5100);
        victim.addPosition(5050);

        Assert.assertEquals(1, victim.build().size());
        assertRegion( victim.build().get(0), 5000, 5100);
    }

    @Test
    public void testOutsideMinGapBehind() {

        final GenomeRegionBuilder victim = new GenomeRegionBuilder(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(5101);

        Assert.assertEquals(2, victim.build().size());
        assertRegion( victim.build().get(0), 5000, 5000);
        assertRegion( victim.build().get(1), 5101, 5101);
    }

    @Test
    public void testJoinMinGapBehind() {

        final GenomeRegionBuilder victim = new GenomeRegionBuilder(CHROM, 100);

        victim.addPosition(5000);
        victim.addPosition(5200);
        victim.addPosition(5100);

        Assert.assertEquals(1, victim.build().size());
        assertRegion( victim.build().get(0), 5000, 5200);
    }

    private void assertRegion(@NotNull final GenomeRegion region, long expectedStart, long expectedEnd) {
        Assert.assertEquals(expectedStart, region.start());
        Assert.assertEquals(expectedEnd, region.end());
    }

}
