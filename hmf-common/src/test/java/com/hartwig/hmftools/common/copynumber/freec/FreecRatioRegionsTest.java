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
        FreecRatio ratio = create("1", 1001, 0.5, 1);
        List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(Lists.newArrayList(ratio));
        assertEquals(1, regions.size());
        assertRegion(regions.get(0), 1001, 2000);
    }

    @Test
    public void testMulipleSimilarRegions() {
        FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(Lists.newArrayList(ratio1, ratio2));
        assertEquals(1, regions.size());
        assertRegion(regions.get(0), 1001, 3000);
    }

    @Test
    public void testMultipleRegionsDifferentBAF() {
        FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        FreecRatio ratio3 = create("1", 3001, 0.7, 1);
        List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(2, regions.size());
        assertRegion(regions.get(0), 1001, 2000);
        assertRegion(regions.get(1), 3001, 4000);
    }

    @Test
    public void testMultipleRegionsMedianRatio() {
        FreecRatio ratio1 = create("1", 1001, 0.5, 1);
        FreecRatio ratio2 = create("1", 2001, 0.5, 1);
        FreecRatio ratio3 = create("1", 3001, 0.5, 1.1);
        List<GenomeRegion> regions = FreecRatioRegions.createRegionsFromRatios(Lists.newArrayList(ratio1, ratio2, ratio3));
        assertEquals(2, regions.size());
        assertRegion(regions.get(0), 1001, 2000);
        assertRegion(regions.get(1), 3001, 4000);
    }

    private void assertRegion(GenomeRegion victim, long expectedStart, long expectedEnd) {
        assertEquals(expectedStart, victim.start());
        assertEquals(expectedEnd, victim.end());
    }


    private static FreecRatio create(@NotNull  String chromosome, long position, double estimatedBaf, double medianRatio) {
        return ImmutableFreecRatio.builder()
                .chromosome(chromosome)
                .position(position)
                .estimatedBAF(estimatedBaf)
                .ratio(medianRatio)
                .medianRatio(medianRatio)
                .build();
    }
}
