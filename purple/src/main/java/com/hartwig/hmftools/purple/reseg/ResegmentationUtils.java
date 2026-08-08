package com.hartwig.hmftools.purple.reseg;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_MAX_PEAKS;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_MAX_SEGMENTATION_PENALTY_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_GAP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_PROPORTION;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_PEAK_MIN_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_STEP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_ROLLING_AVG_WINDOW;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TROUGH_MIN_DIFF;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TROUGH_MIN_GAP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TROUGH_MIN_RATIO;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public class ResegmentationUtils
{
    private final List<RegionData> mOriginalObservedRegions;

    private double mSegmentationPenalty;

    public ResegmentationUtils(final List<ObservedRegion> observedRegions)
    {
        mOriginalObservedRegions = observedRegions.stream().map(x -> new RegionData(x)).collect(Collectors.toList());
        mSegmentationPenalty = 0;
    }

    private class RegionData
    {
        public ObservedRegion OriginalData;

        public double RatioChange;

        public RegionData(final ObservedRegion observedRegion)
        {
            OriginalData = observedRegion;
            RatioChange = 0;
        }

        public String toString() { return OriginalData.toString(); }
    }

    public void processSegments()
    {
        findObservedTumorRatioPenalty();
        findObservedTumorRatioPeak();
    }

    private void findObservedTumorRatioPenalty()
    {
        int ratioCount = (int)((RESEG_RATIO_BUCKET_MAX - RESEG_RATIO_BUCKET_MIN) / RESEG_RATIO_BUCKET_STEP + 1);
        List<TumorRatioDataCount> tumorRatioDataCounts = Lists.newArrayListWithCapacity(ratioCount);

        for(double ratio = 0; ratio <= RESEG_RATIO_BUCKET_MAX; ratio += RESEG_RATIO_BUCKET_STEP)
        {
            tumorRatioDataCounts.add(new TumorRatioDataCount(ratio));
        }

        for(int i = 1; i < mOriginalObservedRegions.size(); ++i)
        {
            RegionData prevRegion = mOriginalObservedRegions.get(i - 1);
            RegionData region = mOriginalObservedRegions.get(i);

            if(region.OriginalData.observedTumorRatio() <= 0 || prevRegion.OriginalData.observedTumorRatio() <= 0)
                continue;

            double ratioChange = abs(region.OriginalData.observedTumorRatio() - prevRegion.OriginalData.observedTumorRatio());
            ratioChange = min(round(ratioChange, 2), RESEG_RATIO_BUCKET_MAX);

            int ratioIndex = (int)(ratioChange / RESEG_RATIO_BUCKET_STEP);
            ++tumorRatioDataCounts.get(ratioIndex).Count;
        }

        calcRollingAverages(ratioCount, tumorRatioDataCounts);

        removeCountDuplicates(tumorRatioDataCounts);

        // find peaks and troughs in the ratio counts
        List<Integer> peakIndices = Lists.newArrayList();
        List<Integer> troughIndices = Lists.newArrayList();
        findPeaksAndTroughs(tumorRatioDataCounts, peakIndices, troughIndices);

        List<TroughData> validTroughs = Lists.newArrayList();

        // find significant troughs vs near-by peaks and troughs
        for(Integer troughIndex : troughIndices)
        {
            double troughAvg = tumorRatioDataCounts.get(troughIndex).RollingAverage;

            Integer prevValidPeakIndex = null;
            Integer nextValidPeakIndex = null;

            // search in each direction
            for(int i = peakIndices.size() - 1; i >= 0; --i)
            {
                if(peakIndices.get(i) > troughIndex)
                    continue;

                int peakIndex = peakIndices.get(i);
                double peakAvg = tumorRatioDataCounts.get(peakIndex).RollingAverage;

                if(peakAvg - troughAvg >= RESEG_TROUGH_MIN_DIFF
                && (troughAvg == 0 || peakAvg / troughAvg >= RESEG_TROUGH_MIN_RATIO))
                {
                    prevValidPeakIndex = peakIndex;
                    break;
                }
            }

            for(int i = 0; i < peakIndices.size(); ++i)
            {
                if(peakIndices.get(i) < troughIndex)
                    continue;

                int peakIndex = peakIndices.get(i);
                double peakAvg = tumorRatioDataCounts.get(peakIndex).RollingAverage;

                if(peakAvg - troughAvg >= RESEG_TROUGH_MIN_DIFF
                && (troughAvg == 0 || peakAvg / troughAvg >= RESEG_TROUGH_MIN_RATIO))
                {
                    nextValidPeakIndex = peakIndex;
                    break;
                }
            }

            if(prevValidPeakIndex != null && nextValidPeakIndex != null)
            {
                validTroughs.add(new TroughData(
                        tumorRatioDataCounts.get(troughIndex).Ratio, troughAvg,
                        tumorRatioDataCounts.get(prevValidPeakIndex).Ratio, tumorRatioDataCounts.get(nextValidPeakIndex).Ratio));
            }
        }

        // consolidate similar troughs
        while(true)
        {
            boolean removedTrough = false;

            int index = 0;
            while(index < validTroughs.size() - 1)
            {
                TroughData lower = validTroughs.get(index);
                TroughData upper = validTroughs.get(index + 1);

                if(upper.Ratio - lower.Ratio >= RESEG_TROUGH_MIN_GAP)
                {
                    ++index;
                }
                else
                {
                    removedTrough = true;

                    if(upper.Average > lower.Average)
                    {
                        validTroughs.remove(index + 1);
                    }
                    else
                    {
                        validTroughs.remove(index);
                    }
                }
            }

            if(!removedTrough)
                break;
        }

        for(TroughData troughData : validTroughs)
        {
            PPL_LOGGER.debug("segmentation trough: {}", troughData);
        }

        // take the lowest peak
        double lowestTroughRatio = validTroughs.stream().mapToDouble(x -> x.Ratio).min().orElse(0);

        lowestTroughRatio = min(lowestTroughRatio, RESEG_MAX_SEGMENTATION_PENALTY_RATIO);

        mSegmentationPenalty = round(pow(lowestTroughRatio, 2) * 2, 4);

        PPL_LOGGER.debug(format("segmentation penalty(%.4f) from lowestTroughRatio(%.f2)",
                mSegmentationPenalty, lowestTroughRatio));
    }

    private class TumorRatioDataCount
    {
        public final double Ratio;
        public int Count;
        public double RollingAverage;

        public TumorRatioDataCount(final double ratio)
        {
            Ratio = ratio;
            Count = 0;
        }

        public String toString() { return format("ratio(%.2f) count(%d) ravg(%.1f)", Ratio, Count, RollingAverage); }
    }

    private class TroughData
    {
        public final double Ratio;
        public final double Average;
        public final double PrevPeakRatio;
        public final double NextPeakRatio;

        public TroughData(final double ratio, final double average, final double prevPeakRatio, final double nextPeakRatio)
        {
            Ratio = ratio;
            Average = average;
            PrevPeakRatio = prevPeakRatio;
            NextPeakRatio = nextPeakRatio;
        }
        public String toString()
        {
            return format("trough(%.2f=%.1f) peaks(prev=%.2f, next=%.2f", Ratio, Average, PrevPeakRatio, NextPeakRatio);
        }
    }

    private class PeakData
    {
        public final double Ratio;
        public final double Average;
        public final double PrevPeakRatio;
        public final double PrevPeakAverage;
        public final double NextPeakRatio;
        public final double NextPeakAverage;

        public PeakData(
                final double ratio, final double average, final double prevPeakRatio, final double prevPeakAverage,
                final double nextPeakRatio, final double nextPeakAverage)
        {
            Ratio = ratio;
            Average = average;
            PrevPeakRatio = prevPeakRatio;
            PrevPeakAverage = prevPeakAverage;
            NextPeakRatio = nextPeakRatio;
            NextPeakAverage = nextPeakAverage;
        }

        public String toString()
        {
            return format("peak(%.2f=%.1f) prev(%.2f=%.1f) next(%.2f=%.1f",
                    Ratio, Average, PrevPeakRatio, PrevPeakAverage, NextPeakRatio, NextPeakRatio);
        }
    }


    private static boolean validPeakRegion(final RegionData region)
    {
        return region.OriginalData.germlineStatus() == GermlineStatus.DIPLOID && region.OriginalData.bafCount() > 0;
    }

    private static void calcRollingAverages(int ratioCount, final List<TumorRatioDataCount> tumorRatioDataCounts)
    {
        int halfRolling = RESEG_ROLLING_AVG_WINDOW / 2;
        for(int i = halfRolling; i < ratioCount - halfRolling; ++i)
        {
            double countTotal = 0;
            for(int k = i - halfRolling; k <= i + halfRolling; ++k)
            {
                countTotal += tumorRatioDataCounts.get(k).Count;
            }

            tumorRatioDataCounts.get(i).RollingAverage = countTotal / RESEG_ROLLING_AVG_WINDOW;
        }
    }

    private static void removeCountDuplicates(final List<TumorRatioDataCount> tumorRatioDataCounts)
    {
        // remove duplicates if any
        int index = 0;
        while(index < tumorRatioDataCounts.size())
        {
            double currentValue = tumorRatioDataCounts.get(index).RollingAverage;

            int nextIndex = index + 1;

            while(nextIndex < tumorRatioDataCounts.size())
            {
                if(tumorRatioDataCounts.get(nextIndex).RollingAverage != currentValue)
                    break;

                tumorRatioDataCounts.remove(nextIndex);
            }

            ++index;
        }
    }

    private static void findPeaksAndTroughs(
            final List<TumorRatioDataCount> tumorRatioDataCounts, final List<Integer> peakIndices, final List<Integer> troughIndices)
    {
        int lastIndex = tumorRatioDataCounts.size()  - 1;

        for(int i = 0; i < tumorRatioDataCounts.size(); ++i)
        {
            double currentValue = tumorRatioDataCounts.get(i).RollingAverage;
            Double prevValue = i > 0 ? tumorRatioDataCounts.get(i - 1).RollingAverage : null;
            Double nextValue = i < lastIndex ? tumorRatioDataCounts.get(i + 1).RollingAverage : null;

            if(prevValue != null && nextValue != null)
            {
                if(currentValue > prevValue && currentValue > nextValue)
                    peakIndices.add(i);
                else if(currentValue < prevValue && currentValue < nextValue)
                    troughIndices.add(i);
            }
            else if(i == 0)
            {
                if(currentValue < nextValue)
                    troughIndices.add(0);
                else if(currentValue > nextValue)
                    peakIndices.add(0);
            }
        }
    }

    private void findObservedTumorRatioPeak()
    {
        int ratioCount = (int)((RESEG_RATIO_BUCKET_MAX - RESEG_RATIO_BUCKET_MIN) / RESEG_RATIO_BUCKET_STEP + 1);
        List<TumorRatioDataCount> tumorRatioDataCounts = Lists.newArrayListWithCapacity(ratioCount);

        for(double ratio = 0; ratio <= RESEG_RATIO_BUCKET_MAX; ratio += RESEG_RATIO_BUCKET_STEP)
        {
            tumorRatioDataCounts.add(new TumorRatioDataCount(ratio));
        }

        // limit regions to diploid and observed BAF counts
        for(int i = 1; i < mOriginalObservedRegions.size(); ++i)
        {
            RegionData region = mOriginalObservedRegions.get(i);

            if(!validPeakRegion(region))
                continue;

            // use same bounds as for tumor ratio diff
            double tumorRatio = max(min(region.OriginalData.observedTumorRatio(), RESEG_RATIO_BUCKET_MAX), 0);

            int ratioIndex = (int)(tumorRatio / RESEG_RATIO_BUCKET_STEP);
            ++tumorRatioDataCounts.get(ratioIndex).Count;
        }

        calcRollingAverages(ratioCount, tumorRatioDataCounts);

        removeCountDuplicates(tumorRatioDataCounts);

        // find peaks and troughs in the ratio counts
        List<Integer> peakIndices = Lists.newArrayList();
        findPeaksAndTroughs(tumorRatioDataCounts, peakIndices, Lists.newArrayList());

        double totalCounts = tumorRatioDataCounts.stream().mapToDouble(x -> x.RollingAverage).sum();

        List<PeakData> validPeaks = Lists.newArrayList();

        // find significant peaks vs near-by peaks
        for(int i = 0; i < peakIndices.size(); ++i)
        {
            Integer peakIndex = peakIndices.get(i);
            double peakAvg = tumorRatioDataCounts.get(peakIndex).RollingAverage;

            if(peakAvg / totalCounts < RESEG_PEAK_MIN_PROPORTION)
                continue;

            Integer prevValidPeakIndex = null;
            Integer nextValidPeakIndex = null;

            // search in each direction
            for(int ud = 0; ud <= 1; ++ud)
            {
                boolean searchUp = (i == 0);
                int k = searchUp ? i + 1 : i - 1;

                for(; k >= 0 && k < peakIndices.size();)
                {
                    int otherPeakIndex = peakIndices.get(k);
                    double otherPeakAvg = tumorRatioDataCounts.get(otherPeakIndex).RollingAverage;

                    if(otherPeakAvg == 0 || peakAvg / otherPeakAvg >= RESEG_PEAK_MIN_RATIO)
                    {
                        if(ud == 0)
                            prevValidPeakIndex = otherPeakIndex;
                        else
                            nextValidPeakIndex = otherPeakIndex;

                        break;
                    }

                    k += searchUp ? 1 : -1;
                }
            }

            if(prevValidPeakIndex != null && nextValidPeakIndex != null)
            {
                PeakData peakData = new PeakData(
                        tumorRatioDataCounts.get(peakIndex).Ratio, peakAvg,
                        tumorRatioDataCounts.get(prevValidPeakIndex).Ratio, tumorRatioDataCounts.get(prevValidPeakIndex).Ratio,
                        tumorRatioDataCounts.get(nextValidPeakIndex).Ratio, tumorRatioDataCounts.get(nextValidPeakIndex).Ratio);

                validPeaks.add(peakData);
            }
        }

        // consolidate similar peak
        while(true)
        {
            boolean removedPeak = false;

            int index = 0;
            while(index < validPeaks.size() - 1)
            {
                PeakData lower = validPeaks.get(index);
                PeakData upper = validPeaks.get(index + 1);

                // [bins[peak_idx] = 0
                // bins[prev_trough] = 2
                // bins[next_trough] = 3
                // peaks.append([bins[peak_idx], values[peak_idx], bins[prev_trough], bins[next_trough], values[prev_trough], values[next_trough]])

                // if right_peak[0] - left_peak[3] >= MIN_PEAK_GAP and right_peak[2] - left_peak[0] >= MIN_PEAK_GAP:

                if(upper.Ratio - lower.NextPeakRatio >= RESEG_PEAK_MIN_GAP
                && upper.PrevPeakRatio - lower.Ratio >= RESEG_PEAK_MIN_GAP)
                {
                    ++index;
                }
                else
                {
                    removedPeak = true;

                    // keep to higher
                    if(upper.Average > lower.Average)
                    {
                        validPeaks.remove(index);
                    }
                    else
                    {
                        validPeaks.remove(index + 1);
                    }
                }
            }

            if(!removedPeak)
                break;
        }

        if(validPeaks.size() > RESEG_MAX_PEAKS)
        {
            List<PeakData> sortedPeaks = Lists.newArrayList(validPeaks);
            Collections.sort(sortedPeaks, Comparator.comparingDouble(x -> -x.Average));
            double minAverage = sortedPeaks.get(RESEG_MAX_PEAKS - 1).Average;

            validPeaks = validPeaks.stream().filter(x -> x.Average >= minAverage).collect(Collectors.toList());
        }

        for(PeakData peakData : validPeaks)
        {
            PPL_LOGGER.debug("segmentation peak: {}", peakData);
        }

        // take the peak with the highest count and then blend in adjacent peaks

        if(validPeaks.size() > 1)
        {
            double maxAverage = validPeaks.get(RESEG_MAX_PEAKS - 1).Average;

            for(int i = 0; i < validPeaks.size(); ++i)
            {
                PeakData peak = validPeaks.get(i);

                if(peak.Average != maxAverage)
                    continue;

                PeakData prevPeak = i > 0 ? validPeaks.get(i - 1) : null;
                PeakData nextPeak = i < validPeaks.size() - 1 ? validPeaks.get(i + 1) : null;

                // To find the left bound, look for the closest peak to the left of the primary peak (by bucket).
                // The left bound is a 2:1 weighted average of the primary peak’s left trough to the nearby peak’s right trough, if the nearby peak exists

            }
        }

        double lowestTroughRatio = validPeaks.stream().mapToDouble(x -> x.Ratio).min().orElse(0);

        lowestTroughRatio = min(lowestTroughRatio, RESEG_MAX_SEGMENTATION_PENALTY_RATIO);

        mSegmentationPenalty = round(pow(lowestTroughRatio, 2) * 2, 4);

        PPL_LOGGER.debug(format("segmentation penalty(%.4f) from lowestTroughRatio(%.f2)",
                mSegmentationPenalty, lowestTroughRatio));
    }

}
