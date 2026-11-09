package com.hartwig.hmftools.purple.reseg;

import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_MAX_PEAKS;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_GAP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_PROPORTION;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MIN;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.reseg.RatioBucketSeries.Bucket;

public final class RatioPeakAnalyser
{
    // finds the per-sample observedTumorRatio peaks (weighted by bafCount, restricted to diploid
    // segments with a positive bafCount) used to scope the GC normalisation
    private RatioPeakAnalyser() {}

    public static RatioPeakResult findRatioPeak(final List<ObservedRegion> segments)
    {
        RatioBucketSeries rawTumorRatios = new RatioBucketSeries();

        for(ObservedRegion segment : segments)
        {
            if(segment.germlineStatus() != GermlineStatus.DIPLOID || segment.bafCount() <= 0)
                continue;

            double ratio = round(segment.observedTumorRatio(), 2);
            rawTumorRatios.addCount(ratio, segment.bafCount());
        }

        List<Bucket> tumorRatios = rawTumorRatios.buildSmoothedSeries();

        if(tumorRatios.size() < 2)
            return null;

        List<Integer> peakIndices = PeakTroughFinder.findLocalPeakOrTroughIndices(tumorRatios, true);

        double totalWeight = tumorRatios.stream().mapToDouble(Bucket::value).sum();

        List<Integer> allIndices = IntStream.range(0, tumorRatios.size()).boxed().collect(Collectors.toList());

        List<PeakTroughData> validPeaks = Lists.newArrayList();

        for(int peakIndex : peakIndices)
        {
            double peakValue = tumorRatios.get(peakIndex).value();

            if(totalWeight <= 0 || peakValue / totalWeight < RESEG_PEAK_MIN_PROPORTION)
                continue;

            OptionalInt leftTrough = PeakTroughFinder.findNearestBelowRequired(
                    allIndices, peakIndex, otherIndex -> isValidPeakSupport(peakValue, tumorRatios.get(otherIndex).value()));

            OptionalInt rightTrough = PeakTroughFinder.findNearestAboveRequired(
                    allIndices, peakIndex, otherIndex -> isValidPeakSupport(peakValue, tumorRatios.get(otherIndex).value()));

            if(leftTrough.isEmpty() || rightTrough.isEmpty())
                continue;

            validPeaks.add(new PeakTroughData(
                    tumorRatios.get(peakIndex).level(), peakValue,
                    tumorRatios.get(leftTrough.getAsInt()).level(), tumorRatios.get(rightTrough.getAsInt()).level()));
        }

        List<PeakTroughData> consolidatedPeaks = PeakTroughFinder.consolidateResults(
                validPeaks,
                (left, right) -> right.Level - left.SupportAboveLevel >= RESEG_PEAK_MIN_GAP
                        && right.SupportBelowLevel - left.Level >= RESEG_PEAK_MIN_GAP,
                x -> x.Value, true);

        List<PeakTroughData> topPeaks = consolidatedPeaks.stream()
                .sorted(Comparator.comparingDouble((PeakTroughData x) -> x.Value).reversed())
                .limit(RESEG_MAX_PEAKS)
                .collect(Collectors.toList());

        if(topPeaks.size() < 2)
            return null;

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

        RatioPeakResult peak = new RatioPeakResult(primary, leftBound, rightBound);
        PPL_LOGGER.trace("primary ratio peak: {}", peak);

        return peak;
    }

    private static boolean isValidPeakSupport(double peakValue, double otherValue)
    {
        return otherValue == 0 || peakValue / otherValue >= RESEG_PEAK_MIN_RATIO;
    }
}
