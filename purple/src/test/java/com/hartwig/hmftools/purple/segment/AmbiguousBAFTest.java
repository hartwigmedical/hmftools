package com.hartwig.hmftools.purple.segment;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.purple.FittingTestUtils.createObservedRegion;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.fitting.AmbiguousBaf;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.junit.Test;

public class AmbiguousBAFTest
{
    private static final double EPSILON = 0.001;

    @Test
    public void testDepthBasedBAF()
    {
        assertEquals(0.581, AmbiguousBaf.depthBasedBAF(30), EPSILON);
        assertEquals(0.5564, AmbiguousBaf.depthBasedBAF(60), EPSILON);
        assertEquals(0.546, AmbiguousBaf.depthBasedBAF(90), EPSILON);
        assertEquals(0.5398, AmbiguousBaf.depthBasedBAF(120), EPSILON);
    }

    @Test
    public void testAmbiguousBAF()
    {
        double validTumorRatio = 1;
        double invalidTumorRatio = 3.1;
        int bafCount = 100;

        // test 1: a single peak
        List<ObservedRegion> observedRegions = List.of(
                createBafObservedRegion(0.49, validTumorRatio, bafCount),
                createBafObservedRegion(0.50, invalidTumorRatio, bafCount),
                createBafObservedRegion(0.51, validTumorRatio, bafCount * 2),
                createBafObservedRegion(0.52, validTumorRatio, bafCount * 3),
                createBafObservedRegion(0.53, validTumorRatio, bafCount * 4),
                createBafObservedRegion(0.54, validTumorRatio, bafCount * 6),
                createBafObservedRegion(0.55, validTumorRatio, bafCount * 4),
                createBafObservedRegion(0.56, validTumorRatio, bafCount * 2),
                createBafObservedRegion(0.57, validTumorRatio, bafCount),
                createBafObservedRegion(0.58, validTumorRatio, bafCount));

        int averageTumorDepth = 100;
        double ambiguousBaf = AmbiguousBaf.findAmbiguousBaf(observedRegions, averageTumorDepth);
        double depthBasedBaf = AmbiguousBaf.depthBasedBAF(averageTumorDepth);

        assertEquals(0.548, ambiguousBaf, EPSILON);

        // no valid points, revert to standard BAF
        observedRegions = List.of(
                createBafObservedRegion(0.52, validTumorRatio, bafCount - 1),
                createBafObservedRegion(0.53, invalidTumorRatio, bafCount),
                createBafObservedRegion(0.54, validTumorRatio, bafCount),
                createBafObservedRegion(0.55, validTumorRatio, bafCount - 1),
                createBafObservedRegion(0.56, invalidTumorRatio, bafCount),
                createBafObservedRegion(0.57, validTumorRatio, bafCount - 1),
                createBafObservedRegion(0.58, validTumorRatio, bafCount - 1));

        ambiguousBaf = AmbiguousBaf.findAmbiguousBaf(observedRegions, averageTumorDepth);
        assertEquals(ambiguousBaf, depthBasedBaf, EPSILON);

        // 2 peaks remaining unconsolidated
        observedRegions = Lists.newArrayList(
                createBafObservedRegion(0.51, validTumorRatio, bafCount),
                createBafObservedRegion(0.515, validTumorRatio, bafCount),
                createBafObservedRegion(0.575, validTumorRatio, bafCount),
                createBafObservedRegion(0.58, validTumorRatio, bafCount));

        ObservedRegion firstPeakRegion = createBafObservedRegion(0.53, validTumorRatio, bafCount * 4);
        buildRegionsAroundPeak(observedRegions, firstPeakRegion, 10);

        ObservedRegion secondPeakRegion = createBafObservedRegion(0.55, validTumorRatio, bafCount * 5);
        buildRegionsAroundPeak(observedRegions, secondPeakRegion, 10);

        Collections.sort(observedRegions, Comparator.comparing(x -> x.observedBAF()));

        ambiguousBaf = AmbiguousBaf.findAmbiguousBaf(observedRegions, averageTumorDepth);

        assertEquals(0.55, ambiguousBaf, EPSILON);

        // revert to using average depth
        depthBasedBaf = AmbiguousBaf.depthBasedBAF(40);

        ambiguousBaf = AmbiguousBaf.findAmbiguousBaf(observedRegions, 40);

        assertEquals(depthBasedBaf, ambiguousBaf, EPSILON);
    }

    private static void buildRegionsAroundPeak(final List<ObservedRegion> observedRegions, final ObservedRegion region, int bafCountDiff)
    {
        int peakWidth = 10;
        int halfWidth = peakWidth / 2;
        double bafDiff = 0.001;

        for(int i = 0; i <= peakWidth; ++i)
        {
            double baf = region.observedBAF() - (halfWidth - i) * bafDiff;
            int depthDiff = abs(halfWidth - i) * bafCountDiff;
            int depth = region.bafCount() - depthDiff;
            ObservedRegion observedRegion = createBafObservedRegion(baf, 2, depth);
            observedRegions.add(observedRegion);
        }

    }

    private static ObservedRegion createBafObservedRegion(double observedBaf, double observedTumorRatio, int bafCount)
    {
        double tumorCN = 1;
        int posStart = 1001;
        int posEnd = 200;

        ObservedRegion observedRegion = createObservedRegion(
                CHR_1, posStart, posEnd, observedBaf, observedTumorRatio, GermlineStatus.DIPLOID, tumorCN);

        observedRegion.setBafCount(bafCount);
        return observedRegion;
    }
}
