package com.hartwig.hmftools.purple.segment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.OptionalInt;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.reseg.PeakTroughData;
import com.hartwig.hmftools.purple.reseg.PeakTroughFinder;
import com.hartwig.hmftools.purple.reseg.PenaltyCalculator;
import com.hartwig.hmftools.purple.reseg.RatioBucketSeries;
import com.hartwig.hmftools.purple.reseg.RatioBucketSeries.Bucket;
import com.hartwig.hmftools.purple.reseg.RatioPeakAnalyser;
import com.hartwig.hmftools.purple.reseg.RatioPeakResult;

import org.junit.Ignore;
import org.junit.Test;

public class ResegmentationTest
{
    // ratio binning and penalties
    @Test
    public void testRollingAverageOnLinearRamp()
    {
        RatioBucketSeries series = new RatioBucketSeries();

        for(int i = 0; i <= 20; i++)
        {
            series.addCount(i * 0.01, i);
        }

        List<Bucket> smoothed = series.buildSmoothedSeries();

        Bucket atTen = smoothed.stream().filter(b -> b.level() == 0.10).findFirst().orElseThrow();
        assertEquals(10.0, atTen.value(), 1e-9);

        Bucket atFive = smoothed.stream().filter(b -> b.level() == 0.05).findFirst().orElseThrow();
        assertEquals(5.0, atFive.value(), 1e-9);
    }

    @Test
    public void testPlateauCollapse()
    {
        RatioBucketSeries series = new RatioBucketSeries();

        for(int i = 50; i <= 60; i++)
        {
            series.addCount(i * 0.01, 10);
        }

        List<Bucket> smoothed = series.buildSmoothedSeries();

        long countOfTens = smoothed.stream().filter(b -> b.value() == 10.0).count();
        assertEquals(1, countOfTens);

        Bucket theTen = smoothed.stream().filter(b -> b.value() == 10.0).findFirst().orElseThrow();
        assertEquals(0.52, theTen.level(), 1e-9);

        // ramps on either side of the flat run remain distinct, non-adjacent occurrences of the same
        // value (6.0 at bucket 0.50 and again at 0.60) are NOT collapsed, since collapsing only merges
        // immediately-adjacent equal values
        long countOfSixes = smoothed.stream().filter(b -> b.value() == 6.0).count();
        assertEquals(2, countOfSixes);
    }

    @Test
    public void testNoAdjacentDuplicateValuesRemainAfterCollapse()
    {
        RatioBucketSeries series = new RatioBucketSeries();

        for(int i = 0; i <= 300; i++)
        {
            series.addCount(i * 0.01, 3);
        }

        List<Bucket> smoothed = series.buildSmoothedSeries();

        // an all-constant raw series collapses down to a single bucket
        assertEquals(1, smoothed.size());
        assertTrue(smoothed.get(0).value() == 3.0);
    }

    // peaks and trough testing
    @Ignore
    @Test
    public void testTwoClusterDiploidSegmentsYieldsPrimaryPeakAndBounds()
    {
        List<ObservedRegion> segments = new ArrayList<>();
        int pos = 1;

        // cluster around observedTumorRatio 1.00, weighted by bafCount 200 -> smoothed plateau value 40
        for(int i = 0; i < 50; i++)
        {
            segments.add(bothNoneSegment("1", pos, pos + 999, 200, 1.00));
            pos += 1000;
        }

        // smaller cluster around observedTumorRatio 1.50, bafCount 100 -> smoothed plateau value 20
        for(int i = 0; i < 50; i++)
        {
            segments.add(bothNoneSegment("1", pos, pos + 999, 100, 1.50));
            pos += 1000;
        }

        Optional<RatioPeakResult> result = RatioPeakAnalyser.findRatioPeak(segments);

        assertTrue(result.isPresent());

        // the smoothing window spreads each ratio spike into a 5-wide plateau; the peak level reported
        // is the first (lowest-level) bucket of that plateau, not the raw spike ratio itself
        double primaryLevel = result.get().PrimaryPeak.Level;
        assertEquals(0.98, primaryLevel, 1e-9);
        assertTrue(result.get().LeftBound < primaryLevel);
        assertTrue(result.get().RightBound > primaryLevel);
    }

    @Test
    public void testSingleClusterYieldsNoResult()
    {
        List<ObservedRegion> segments = new ArrayList<>();
        int pos = 1;

        for(int i = 0; i < 50; i++)
        {
            segments.add(bothNoneSegment("1", pos, pos + 999, 200, 1.00));
            pos += 1000;
        }

        Optional<RatioPeakResult> result = RatioPeakAnalyser.findRatioPeak(segments);

        assertFalse(result.isPresent());
    }

    // tumor ratio penalty
    @Test
    public void testEngineeredBimodalDistributionYieldsExpectedTrough()
    {
        List<ObservedRegion> segments = new ArrayList<>();
        int pos = 1;

        // cluster 1: 101 segments alternating ratio 1.00/1.10 -> 100 adjacent diffs of exactly 0.10
        // (bucket index 10, smoothed plateau value 100/5 = 20)
        for(int i = 0; i <= 100; i++)
        {
            double ratio = (i % 2 == 0) ? 1.00 : 1.10;
            segments.add(bothNoneSegment("1", pos, pos + 999, 100, ratio));
            pos += 1000;
        }

        // cluster 2: 101 segments alternating ratio 4.00/4.40 -> 100 adjacent diffs of exactly 0.40
        // (bucket index 40, smoothed plateau value 100/5 = 20). The single cross-cluster diff
        // (1.00 -> 4.00 = 3.00) clamps to the far-away max bucket and doesn't interact with either plateau.
        for(int i = 0; i <= 100; i++)
        {
            double ratio = (i % 2 == 0) ? 4.00 : 4.40;
            segments.add(bothNoneSegment("1", pos, pos + 999, 100, ratio));
            pos += 1000;
        }

        double penalty = PenaltyCalculator.calculatePenalty(segments);

        // expected trough sits at the first bucket of the flat zero run between the two plateaus (level 0.13):
        // penalty = 0.13^2 * 2 = 0.0338
        assertEquals(0.0338, penalty, 1e-9);
    }

    @Test
    public void testDegenerateAllNonPositiveRatiosUsesDefaultPenalty()
    {
        List<ObservedRegion> segments = new ArrayList<>();

        for(int i = 0; i < 5; i++)
        {
            segments.add(bothNoneSegment("1", i * 1000 + 1, i * 1000 + 999, 100, 0.0));
        }

        double penalty = PenaltyCalculator.calculatePenalty(segments);

        // no valid ratio-diff pairs at all (all ratios are 0) -> falls back to the 0.35 default bucket
        assertEquals(0.245, penalty, 1e-9);
    }


    private static List<Bucket> bucketsOf(double... values)
    {
        Bucket[] buckets = new Bucket[values.length];
        for(int i = 0; i < values.length; i++)
        {
            buckets[i] = new Bucket(i * 0.01, values[i]);
        }
        return Arrays.asList(buckets);
    }

    @Test
    public void testFindLocalExtremaIndices()
    {
        List<Bucket> series = bucketsOf(1, 3, 2, 5, 1, 4, 2);

        List<Integer> peaks = PeakTroughFinder.findLocalPeakOrTroughIndices(series, true);
        assertEquals(Arrays.asList(1, 3, 5), peaks);

        List<Integer> troughs = PeakTroughFinder.findLocalPeakOrTroughIndices(series, false);
        assertEquals(Arrays.asList(2, 4), troughs);
    }

    @Test
    public void testFindNearestSatisfyingStopsAtFirstMatchMovingOutward()
    {
        List<Bucket> series = bucketsOf(0, 1, 2, 3, 4, 5, 6);
        List<Integer> candidates = Arrays.asList(1, 3, 5);

        // only candidate 3 satisfies -> nearest-below search from center 4 should find it directly
        OptionalInt belowDirect = PeakTroughFinder.findNearestBelowRequired(candidates, 4, idx -> idx == 3);
        assertTrue(belowDirect.isPresent());
        assertEquals(3, belowDirect.getAsInt());

        // only candidate 1 satisfies (not the nearer candidate 3) -> must walk further out
        OptionalInt belowFurther = PeakTroughFinder.findNearestBelowRequired(candidates, 4, idx -> idx == 1);
        assertTrue(belowFurther.isPresent());
        assertEquals(1, belowFurther.getAsInt());

        // only candidate 5 satisfies, searching above center 2 (must skip nearer candidate 3)
        OptionalInt above = PeakTroughFinder.findNearestAboveRequired(candidates, 2, idx -> idx == 5);
        assertTrue(above.isPresent());
        assertEquals(5, above.getAsInt());

        // nothing satisfies
        OptionalInt none = PeakTroughFinder.findNearestBelowRequired(candidates, 4, idx -> false);
        assertFalse(none.isPresent());
    }

    @Test
    public void testConsolidateTroughsKeepsLowerValueOnClose()
    {
        List<PeakTroughData> troughs = Arrays.asList(
                new PeakTroughData(0.10, 5, 0, 0),
                new PeakTroughData(0.15, 3, 0, 0),
                new PeakTroughData(0.30, 1, 0, 0));

        List<PeakTroughData> result = PeakTroughFinder.consolidateResults(
                troughs, (left, right) -> right.Level - left.Level >= 0.1, x -> x.Value, false);

        assertEquals(2, result.size());
        assertEquals(0.15, result.get(0).Level, 1e-9);
        assertEquals(0.30, result.get(1).Level, 1e-9);
    }

    @Test
    public void testConsolidatePeaksKeepsHigherValueOnClose()
    {
        PeakTroughData a = new PeakTroughData(0.10, 5, 0, 0.20);
        PeakTroughData b = new PeakTroughData(0.15, 8, 0.12, 0.25);

        List<PeakTroughData> result = PeakTroughFinder.consolidateResults(
                Arrays.asList(a, b),
                (left, right) -> right.Level - left.SupportAboveLevel >= 0.1 && right.SupportBelowLevel - left.Level >= 0.1,
                x -> x.Value, true);

        assertEquals(1, result.size());
        assertEquals(0.15, result.get(0).Level, 1e-9);
    }

    @Test
    public void testConsolidateIsIdempotent()
    {
        List<PeakTroughData> troughs = Arrays.asList(
                new PeakTroughData(0.10, 5, 0, 0),
                new PeakTroughData(0.15, 3, 0, 0),
                new PeakTroughData(0.30, 1, 0, 0));

        List<PeakTroughData> once = PeakTroughFinder.consolidateResults(
                troughs, (left, right) -> right.Level - left.Level >= 0.1, x -> x.Value, false);

        List<PeakTroughData> twice = PeakTroughFinder.consolidateResults(
                once, (left, right) -> right.Level - left.Level >= 0.1, x -> x.Value, false);

        assertEquals(once.size(), twice.size());
        for(int i = 0; i < once.size(); i++)
        {
            assertEquals(once.get(i).Level, twice.get(i).Level, 1e-9);
        }
    }

    private static ObservedRegion segment(
            final String chromosome, int start, int end, final SegmentSupport support, final GermlineStatus germlineStatus,
            int bafCount, double observedTumorRatio, double gcContent, int depthWindowCount)
    {
        return new ObservedRegion(
                chromosome, start, end, true, support, bafCount, 0.5, depthWindowCount,
                observedTumorRatio, 1.0, 1.0, germlineStatus, false, gcContent,
                start, start, 0, 0, 0, 0,
                0, 0, 0, 0, 0);
    }

    private static ObservedRegion bothNoneSegment(
            final String chromosome, int start, int end, int bafCount, double observedTumorRatio, double gcContent, int depthWindowCount)
    {
        return segment(
                chromosome, start, end, SegmentSupport.NONE, GermlineStatus.DIPLOID,
                bafCount, observedTumorRatio, gcContent, depthWindowCount);
    }

    private static ObservedRegion bothNoneSegment(final String chromosome, int start, int end, int bafCount, double observedTumorRatio)
    {
        return bothNoneSegment(chromosome, start, end, bafCount, observedTumorRatio, 0.45, 100);
    }
}
