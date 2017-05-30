package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;
import static com.hartwig.hmftools.common.numeric.Doubles.lessThan;
import static com.hartwig.hmftools.common.purple.region.SmoothedRegions.allowedBAFDeviation;
import static com.hartwig.hmftools.common.purple.region.SmoothedRegions.allowedCopyNumberDeviation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedRegion;

import org.junit.Test;

public class SmoothedRegionTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void copyNumberAllowances() {
        assertEquals(1.3, allowedCopyNumberDeviation(0), EPSILON);
        assertEquals(0.8, allowedCopyNumberDeviation(5), EPSILON);
        assertEquals(0.3, allowedCopyNumberDeviation(10), EPSILON);
        assertEquals(0.3, allowedCopyNumberDeviation(50), EPSILON);
    }

    @Test
    public void bafAllowances() {
        assertTrue(lessThan(allowedBAFDeviation(1), 0.44));
        assertTrue(greaterThan(allowedBAFDeviation(7), 0.072));
        assertTrue(greaterThan(allowedBAFDeviation(14), 0.052));
        assertTrue(greaterThan(allowedBAFDeviation(46), 0.03));
    }

    @Test
    public void beforeFirstBroadRegion() {
        final List<ConsolidatedRegion> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2));
        final List<FittedRegion> copyNumbers = Lists.newArrayList(
                createFittedCopyNumber(1, 300, 1),
                createFittedCopyNumber(301, 1000, 2),
                createFittedCopyNumber(1001, 2000, 2));

        final List<ConsolidatedRegion> results = new SmoothedRegions(broadRegions, copyNumbers).getSmoothedRegions();
        assertEquals(2, results.size());

        assertRegion(results.get(0), 1, 300, 1);
        assertRegion(results.get(1), 301, 2000, 2);
    }

    @Test
    public void betweenBroadRegions() {
        final List<ConsolidatedRegion> broadRegions = Lists.newArrayList(
                createRegion(1001, 2000, 2),
                createRegion(3001, 4000, 3));
        final List<FittedRegion> copyNumbers = Lists.newArrayList(
                createFittedCopyNumber(1001, 2000, 2),
                createFittedCopyNumber(2001, 2200, 2),
                createFittedCopyNumber(2201, 2500, 1),
                createFittedCopyNumber(2501, 3000, 3),
                createFittedCopyNumber(3001, 4000, 3));

        final List<ConsolidatedRegion> results = new SmoothedRegions(broadRegions, copyNumbers).getSmoothedRegions();
        assertEquals(3, results.size());

        assertRegion(results.get(0), 1001, 2200, 2);
        assertRegion(results.get(1), 2201, 2500, 1);
        assertRegion(results.get(2), 2501, 4000, 3);
    }

    @Test
    public void afterLastBroadRegion() {
        final List<ConsolidatedRegion> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2));
        final List<FittedRegion> copyNumbers = Lists.newArrayList(
                createFittedCopyNumber(1001, 2000, 2),
                createFittedCopyNumber(2001, 3000, 2),
                createFittedCopyNumber(3001, 5000, 5));

        final List<ConsolidatedRegion> results = new SmoothedRegions(broadRegions, copyNumbers).getSmoothedRegions();
        assertEquals(2, results.size());

        assertRegion(results.get(0), 1001, 3000, 2);
        assertRegion(results.get(1), 3001, 5000, 5);
    }


    private void assertRegion(ConsolidatedRegion victim, long expectedStart, long expectedEnd, double expectedCopyNumber) {
        assertEquals(expectedStart, victim.start());
        assertEquals(expectedEnd, victim.end());
        assertEquals(expectedCopyNumber, victim.averageTumorCopyNumber(), EPSILON);
    }

    private FittedRegion createFittedCopyNumber(long start, long end, double copyNumber) {
        return ConsolidatedRegionBuilderTest.create("1", start, end, 10, 0.5, copyNumber);
    }

    private ConsolidatedRegion createRegion(long start, long end, double copyNumber) {
        return ImmutableConsolidatedRegion.builder()
                .chromosome("1")
                .start(start)
                .end(end)
                .averageTumorCopyNumber(copyNumber)
                .bafCount(10)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5)
                .build();
    }
}
