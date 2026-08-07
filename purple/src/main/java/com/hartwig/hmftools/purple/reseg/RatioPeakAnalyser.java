package com.hartwig.hmftools.purple.reseg;

import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_MAX_PEAKS;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_GAP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_PROPORTION;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MIN;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.reseg.RatioBucketSeries.Bucket;

public final class RatioPeakAnalyser
{
    // finds the per-sample observedTumorRatio peaks (weighted by bafCount, restricted to diploid
    // segments with a positive bafCount) used to scope the GC normalisation
    private RatioPeakAnalyser() {}

    public static Optional<RatioPeakResult> analyse(final List<ObservedRegion> segments)
    {
        RatioBucketSeries bucketSeries = new RatioBucketSeries();

        for(ObservedRegion segment : segments)
        {
            if(segment.germlineStatus() != GermlineStatus.DIPLOID || segment.bafCount() <= 0)
                continue;

            double ratio = round(segment.observedTumorRatio(), 2);
            bucketSeries.addCount(ratio, segment.bafCount());
        }

        List<Bucket> series = bucketSeries.buildSmoothedSeries();

        if(series.size() < 2)
            return Optional.empty();

        List<Integer> peakIndices = PeakTroughFinder.findLocalExtremaIndices(series, true);

        double totalWeight = series.stream().mapToDouble(Bucket::value).sum();

        List<Integer> allIndices = IntStream.range(0, series.size()).boxed().collect(Collectors.toList());

        List<PeakTroughData> validPeaks = new ArrayList<>();

        for(int peakIndex : peakIndices)
        {
            double peakValue = series.get(peakIndex).value();

            if(totalWeight <= 0 || peakValue / totalWeight < RESEG_PEAK_MIN_PROPORTION)
                continue;

            OptionalInt leftTrough = PeakTroughFinder.findNearestSatisfyingBelow(
                    allIndices, peakIndex, otherIndex -> isValidPeakSupport(peakValue, series.get(otherIndex).value()));

            OptionalInt rightTrough = PeakTroughFinder.findNearestSatisfyingAbove(
                    allIndices, peakIndex, otherIndex -> isValidPeakSupport(peakValue, series.get(otherIndex).value()));

            if(leftTrough.isEmpty() || rightTrough.isEmpty())
                continue;

            validPeaks.add(new PeakTroughData(
                    series.get(peakIndex).level(), peakValue,
                    series.get(leftTrough.getAsInt()).level(), series.get(rightTrough.getAsInt()).level()));
        }

        List<PeakTroughData> consolidated = PeakTroughFinder.consolidate(
                validPeaks,
                (left, right) -> right.Level - left.SupportAboveLevel >= RESEG_PEAK_MIN_GAP
                        && right.SupportBelowLevel - left.Level >= RESEG_PEAK_MIN_GAP,
                x -> x.Value, true);

        List<PeakTroughData> topPeaks = consolidated.stream()
                .sorted(Comparator.comparingDouble((PeakTroughData x) -> x.Value).reversed())
                .limit(RESEG_MAX_PEAKS)
                .collect(Collectors.toList());

        if(topPeaks.size() < 2)
            return Optional.empty();

        PeakTroughData primary = topPeaks.get(0);
        List<PeakTroughData> others = topPeaks.subList(1, topPeaks.size());

        double leftBound = others.stream()
                .filter(p -> p.Level < primary.Level)
                .max(Comparator.comparingDouble(p -> p.Level))
                .map(p -> round((p.SupportAboveLevel + 2 * primary.SupportBelowLevel) / 3, 4))
                .orElse(RESEG_RATIO_BUCKET_MIN);

        double rightBound = others.stream()
                .filter(p -> p.Level > primary.Level)
                .min(Comparator.comparingDouble(p -> p.Level))
                .map(p -> round((p.SupportBelowLevel + 2 * primary.SupportAboveLevel) / 3, 4))
                .orElse(RESEG_RATIO_BUCKET_MAX);

        return Optional.of(new RatioPeakResult(primary, leftBound, rightBound));
    }

    private static boolean isValidPeakSupport(double peakValue, double otherValue)
    {
        return otherValue == 0 || peakValue / otherValue >= RESEG_PEAK_MIN_RATIO;
    }
}
