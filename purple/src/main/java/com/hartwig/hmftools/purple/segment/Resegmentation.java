package com.hartwig.hmftools.purple.segment;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_DIFF_BOUND_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_DIFF_BOUND_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_DIFF_STEP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_RATIO_PENALTY;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_ROLL_AVG_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_TROUGH_DIFF;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_TROUGH_GAP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_TUMOR_RATIO_TROUGH_RATIO;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public class Resegmentation
{
    private final List<RegionData> mOriginalObservedRegions;

    private double mSegmentationPenalty;

    public Resegmentation(final List<ObservedRegion> observedRegions)
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
    }

    private void findObservedTumorRatioPenalty()
    {
        int ratioCount = (int)((RESEG_TUMOR_RATIO_DIFF_BOUND_MAX - RESEG_TUMOR_RATIO_DIFF_BOUND_MIN) / RESEG_TUMOR_RATIO_DIFF_STEP + 1);
        List<TumorRatioChangeCount> tumorRatioChangeCounts = Lists.newArrayListWithCapacity(ratioCount);

        for(double ratio = 0; ratio <= RESEG_TUMOR_RATIO_DIFF_BOUND_MAX; ratio += RESEG_TUMOR_RATIO_DIFF_STEP)
        {
            tumorRatioChangeCounts.add(new TumorRatioChangeCount(ratio));
        }

        for(int i = 1; i < mOriginalObservedRegions.size(); ++i)
        {
            RegionData prevRegion = mOriginalObservedRegions.get(i - 1);
            RegionData region = mOriginalObservedRegions.get(i);

            if(region.OriginalData.observedTumorRatio() <= 0 || prevRegion.OriginalData.observedTumorRatio() <= 0)
                continue;

            double ratioChange = abs(region.OriginalData.observedTumorRatio() - prevRegion.OriginalData.observedTumorRatio());
            ratioChange = min(round(ratioChange, 2), RESEG_TUMOR_RATIO_DIFF_BOUND_MAX);

            int ratioIndex = (int)(ratioChange / RESEG_TUMOR_RATIO_DIFF_STEP);
            ++tumorRatioChangeCounts.get(ratioIndex).Count;
        }

        int halfRolling = RESEG_TUMOR_RATIO_ROLL_AVG_COUNT / 2;
        for(int i = halfRolling; i < ratioCount - halfRolling; ++i)
        {
            double countTotal = 0;
            for(int k = i - halfRolling; k <= i + halfRolling; ++k)
            {
                countTotal += tumorRatioChangeCounts.get(k).Count;
            }

            tumorRatioChangeCounts.get(i).RollingAverage = countTotal / RESEG_TUMOR_RATIO_ROLL_AVG_COUNT;
        }

        // remove duplicates if any
        int index = 0;
        while(index < tumorRatioChangeCounts.size())
        {
            double currentValue = tumorRatioChangeCounts.get(index).RollingAverage;

            int nextIndex = index + 1;

            while(nextIndex < tumorRatioChangeCounts.size())
            {
                if(tumorRatioChangeCounts.get(nextIndex).RollingAverage != currentValue)
                    break;

                tumorRatioChangeCounts.remove(nextIndex);
            }

            ++index;
        }

        // find peaks and troughs in the ratio counts
        List<Integer> peakIndices = Lists.newArrayList();
        List<Integer> troughIndices = Lists.newArrayList();
        int lastIndex = tumorRatioChangeCounts.size()  - 1;

        for(int i = 0; i < tumorRatioChangeCounts.size(); ++i)
        {
            double currentValue = tumorRatioChangeCounts.get(i).RollingAverage;
            Double prevValue = i > 0 ? tumorRatioChangeCounts.get(i - 1).RollingAverage : null;
            Double nextValue = i < lastIndex ? tumorRatioChangeCounts.get(i + 1).RollingAverage : null;

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

        List<TroughData> validTroughs = Lists.newArrayList();

        // find significant troughs vs near-by peaks and troughs
        for(Integer troughIndex : troughIndices)
        {
            double troughAvg = tumorRatioChangeCounts.get(troughIndex).RollingAverage;

            Integer prevValidPeakIndex = null;
            Integer nextValidPeakIndex = null;

            // search in each direction
            for(int i = peakIndices.size() - 1; i >= 0; --i)
            {
                if(peakIndices.get(i) > troughIndex)
                    continue;

                int peakIndex = peakIndices.get(i);
                double peakAvg = tumorRatioChangeCounts.get(peakIndex).RollingAverage;

                if(peakAvg - troughAvg >= RESEG_TUMOR_RATIO_TROUGH_DIFF
                && (troughAvg == 0 || peakAvg / troughAvg >= RESEG_TUMOR_RATIO_TROUGH_RATIO))
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
                double peakAvg = tumorRatioChangeCounts.get(peakIndex).RollingAverage;

                if(peakAvg - troughAvg >= RESEG_TUMOR_RATIO_TROUGH_DIFF
                && (troughAvg == 0 || peakAvg / troughAvg >= RESEG_TUMOR_RATIO_TROUGH_RATIO))
                {
                    nextValidPeakIndex = peakIndex;
                    break;
                }
            }

            if(prevValidPeakIndex != null && nextValidPeakIndex != null)
            {
                validTroughs.add(new TroughData(
                        tumorRatioChangeCounts.get(troughIndex).Ratio, troughAvg,
                        tumorRatioChangeCounts.get(prevValidPeakIndex).Ratio, tumorRatioChangeCounts.get(nextValidPeakIndex).Ratio));
            }
        }

        // consolidate similar troughs
        while(true)
        {
            boolean removedTrough = false;

            index = 0;
            while(index < validTroughs.size() - 1)
            {
                TroughData lower = validTroughs.get(index);
                TroughData upper = validTroughs.get(index + 1);

                if(upper.Ratio - lower.Ratio >= RESEG_TUMOR_RATIO_TROUGH_GAP)
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

        lowestTroughRatio = min(lowestTroughRatio, RESEG_TUMOR_RATIO_RATIO_PENALTY);

        mSegmentationPenalty = round(pow(lowestTroughRatio, 2) * 2, 4);

        PPL_LOGGER.debug(format("segmentation penalty(%.4f) from lowestTroughRatio(%.f2)",
                mSegmentationPenalty, lowestTroughRatio));
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

    private class TumorRatioChangeCount
    {
        public final double Ratio;
        public int Count;
        public double RollingAverage;

        public TumorRatioChangeCount(final double ratio)
        {
            Ratio = ratio;
            Count = 0;
        }

        public String toString() { return format("ratio(%.2f) count(%d) ravg(%.1f)", Ratio, Count, RollingAverage); }
    }
}
