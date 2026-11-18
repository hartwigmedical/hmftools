package com.hartwig.hmftools.purple.reseg;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_MAX_SEGMENTATION_PENALTY_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TROUGH_MIN_DIFF;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TROUGH_MIN_GAP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TROUGH_MIN_RATIO;

import java.util.List;
import java.util.OptionalInt;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.reseg.RatioBucketSeries.Bucket;

public final class PenaltyCalculator
{
    // calc segmentation penalty from the observedTumorRatio diffs, looking at peaks and troughs
    private PenaltyCalculator() {}

    private static double calcPenality(double ratio)
    {
        return round(Math.pow(ratio, 2) * 2, 4);
    }

    private static final double DEFAULT_PENALTY = calcPenality(RESEG_MAX_SEGMENTATION_PENALTY_RATIO);

    public static double calculatePenalty(final List<ObservedRegion> segmentsInOrder)
    {
        RatioBucketSeries rawTumorRatioDiffs = new RatioBucketSeries();

        for(int i = 1; i < segmentsInOrder.size(); i++)
        {
            ObservedRegion prev = segmentsInOrder.get(i - 1);
            ObservedRegion curr = segmentsInOrder.get(i);

            if(curr.observedTumorRatio() <= 0 || prev.observedTumorRatio() <= 0)
                continue;

            double ratioChange = Math.abs(curr.observedTumorRatio() - prev.observedTumorRatio());
            ratioChange = min(round(ratioChange, 2), RESEG_RATIO_BUCKET_MAX);

            rawTumorRatioDiffs.addCount(ratioChange, 1);
        }

        List<Bucket> tumorRatioDiffs = rawTumorRatioDiffs.buildSmoothedSeries();

        if(tumorRatioDiffs.size() < 2)
            return DEFAULT_PENALTY;

        List<Integer> peakIndices = PeakTroughFinder.findLocalPeakOrTroughIndices(tumorRatioDiffs, true);
        List<Integer> troughIndices = PeakTroughFinder.findLocalPeakOrTroughIndices(tumorRatioDiffs, false);

        // add the first bucket as a peak or trough
        if(tumorRatioDiffs.get(0).value() < tumorRatioDiffs.get(1).value())
            troughIndices.add(0, 0);
        else if(tumorRatioDiffs.get(0).value() > tumorRatioDiffs.get(1).value())
            peakIndices.add(0, 0);

        List<PeakTroughData> validTroughs = Lists.newArrayList();

        for(int troughIndex : troughIndices)
        {
            double troughValue = tumorRatioDiffs.get(troughIndex).value();

            OptionalInt prevPeak = PeakTroughFinder.findNearestBelowRequired(
                    peakIndices, troughIndex, peakIndex -> isValidTroughSupport(troughValue, tumorRatioDiffs.get(peakIndex).value()));

            OptionalInt nextPeak = PeakTroughFinder.findNearestAboveRequired(
                    peakIndices, troughIndex, peakIndex -> isValidTroughSupport(troughValue, tumorRatioDiffs.get(peakIndex).value()));

            if(prevPeak.isEmpty() || nextPeak.isEmpty())
                continue;

            validTroughs.add(new PeakTroughData(
                    tumorRatioDiffs.get(troughIndex).level(), troughValue,
                    tumorRatioDiffs.get(prevPeak.getAsInt()).level(), tumorRatioDiffs.get(nextPeak.getAsInt()).level()));
        }

        List<PeakTroughData> consolidatedTroughs = PeakTroughFinder.consolidateResults(
                validTroughs, (left, right) -> right.Level - left.Level >= RESEG_TROUGH_MIN_GAP, x -> x.Value, false);

        double lowestTroughLevel = consolidatedTroughs.stream().mapToDouble(x -> x.Level).min().orElse(RESEG_MAX_SEGMENTATION_PENALTY_RATIO);
        double penaltyBucket = min(lowestTroughLevel, RESEG_MAX_SEGMENTATION_PENALTY_RATIO);

        return calcPenality(penaltyBucket);
    }

    private static boolean isValidTroughSupport(double troughValue, double peakValue)
    {
        if(peakValue - troughValue < RESEG_TROUGH_MIN_DIFF)
            return false;

        return troughValue == 0 || peakValue / troughValue >= RESEG_TROUGH_MIN_RATIO;
    }
}
