package com.hartwig.hmftools.purple.reseg;

import static java.lang.Double.NaN;
import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Doubles.median;
import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_GC_PERCENTILES;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_GC_SPARSE_MAX_PROPORTION;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_GC_SPARSE_MIN_JUMP;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.Nullable;

public final class GcNormaliser
{
    // per-sample bucketed GC-content normalisation of observedTumorRatio, restricted to diploid
    // segments with a positive bafCount. If a primary ratio peak was found, also restrict to with its bounds.
    private static final int MAX_SUPPRESSION_PASSES = 20;

    public static List<ObservedRegion> normaliseRegions(final List<ObservedRegion> allSegments, @Nullable final RatioPeakResult peakResult)
    {
        List<GcBin> allBins = GcBin.createdBins();

        // the set of bins actually touched by any segment in the full sample, used to resolve the
        // out-of-range LOW/HIGH clamp bins to the observed extremes rather than the theoretical ones
        Set<GcBin> observedBins = Sets.newHashSet();

        GcBin lowerGcBin = null;
        GcBin highGcBin = null;

        for(ObservedRegion segment : allSegments)
        {
            GcBin gcBin = findBin(segment.gcContent(), allBins);
            if(gcBin != null)
            {
                observedBins.add(gcBin);

                if(lowerGcBin == null || gcBin.Lower < lowerGcBin.Lower)
                    lowerGcBin = gcBin;

                if(highGcBin == null || gcBin.Upper > highGcBin.Upper)
                    highGcBin = gcBin;

                if(observedBins.size() == allBins.size())
                    break;
            }
        }

        if(lowerGcBin == null)
            lowerGcBin = allBins.get(0);

        if(highGcBin == null)
            highGcBin = allBins.get(allBins.size() - 1);

        List<ObservedRegion> validIputSegments = allSegments.stream()
                .filter(x -> isValidRegion(x, peakResult))
                .collect(Collectors.toList());

        Map<GcBin,Double> gcAdjustments = calculateGcAdjustments(validIputSegments, allBins);

        List<ObservedRegion> gcAdjustedRegions = Lists.newArrayListWithCapacity(allSegments.size());

        for(ObservedRegion segment : allSegments)
        {
            GcBin bin = findBin(segment.gcContent(), allBins);

            if(bin == null)
                bin = segment.gcContent() <= lowerGcBin.Lower ? lowerGcBin : highGcBin;

            double adjustment = gcAdjustments.getOrDefault(bin, 1.0);

            if(Double.isNaN(adjustment) || adjustment <= 0)
                adjustment = 1.0;

            // double newRatio = round(segment.observedTumorRatio() / adjustment, 4);
            double newRatio = segment.observedTumorRatio() / adjustment; // skip rounding

            ObservedRegion newRegion = ObservedRegion.fromOther(segment);
            newRegion.setObservedTumorRatio(newRatio);
            gcAdjustedRegions.add(newRegion);
        }

        return gcAdjustedRegions;
    }

    private static boolean isValidRegion(final ObservedRegion region, @Nullable final RatioPeakResult peakResult)
    {
        if(region.bafCount() <= 0)
            return false;

        if(region.germlineStatus() != GermlineStatus.DIPLOID)
            return false;

        if(peakResult != null)
        {
            if(region.observedTumorRatio() < peakResult.LeftBound || region.observedTumorRatio() > peakResult.RightBound)
                return false;
        }

        return true;
    }

    private static Map<GcBin,Double> calculateGcAdjustments(final List<ObservedRegion> regions, final List<GcBin> allBins)
    {
        double totalWeight = regions.stream().mapToDouble(s -> s.bafCount()).sum();

        double[] binProportion = new double[allBins.size()];

        List<List<ObservedRegion>> regionsByGcBin = Lists.newArrayListWithCapacity(allBins.size());

        for(int i = 0; i < allBins.size(); i++)
        {
            GcBin bin = allBins.get(i);

            List<ObservedRegion> binRegions = Lists.newArrayList();
            regionsByGcBin.add(binRegions);

            double binWeight = 0;

            for(ObservedRegion region : regions)
            {
                int bafCount = region.bafCount();
                if(bin.contains(region.gcContent()))
                {
                    binWeight += bafCount;
                    binRegions.add(region);
                }
            }

            binProportion[i] = totalWeight > 0 ? binWeight / totalWeight : 0;
        }

        List<RatioBafCount> allSegmentRatioCounts = buildRatioBafCounts(regions);

        double[][] grid = new double[allBins.size()][RESEG_GC_PERCENTILES.length];

        for(int p = 0; p < RESEG_GC_PERCENTILES.length; p++)
        {
            double percentile = RESEG_GC_PERCENTILES[p];

            for(int i = 0; i < allBins.size(); i++)
            {
                List<ObservedRegion> gcSegments = regionsByGcBin.get(i);
                List<RatioBafCount> gcSegmentRatioCounts = buildRatioBafCounts(gcSegments);

                Double allSegsPercRatio = weightRatioAcrossPercentile(allSegmentRatioCounts, percentile);
                Double gcSegsPercRatio = weightRatioAcrossPercentile(gcSegmentRatioCounts, percentile);

                if(allSegsPercRatio != null && gcSegsPercRatio != null && allSegsPercRatio != null)
                    grid[i][p] = gcSegsPercRatio / allSegsPercRatio;
                else
                    grid[i][p] = NaN;
            }
        }

        suppressSparseValues(grid, binProportion);

        Map<GcBin,Double> adjustments = new LinkedHashMap<>();

        for(int i = 0; i < allBins.size(); i++)
        {
            List<Double> validEntries = Lists.newArrayList();

            for(double v : grid[i])
            {
                if(!Double.isNaN(v))
                    validEntries.add(v);
            }

            double medianAdjustment = validEntries.isEmpty() ? NaN : median(validEntries); // no rounding to 4dp as in protoype

            // if(!Double.isNaN(medianAdjustment))
            //    medianAdjustment = round(medianAdjustment, 4);

            adjustments.put(allBins.get(i), medianAdjustment);
        }

        return adjustments;
    }

    private static void suppressSparseValues(final double[][] grid, final double[] binProportion)
    {
        // mark entries that jump too much from the preceding bin and are backed by too little data as unreliable, then forward-fill
        // them from the nearest preceding valid bin, repeating until a full pass finds nothing left to suppress
        int binCount = grid.length;
        int percentileCount = grid[0].length;

        boolean anyNaN = true;
        int pass = 0;

        while(anyNaN && pass < MAX_SUPPRESSION_PASSES)
        {
            pass++;

            for(int p = 0; p < percentileCount; p++)
            {
                for(int i = 1; i < binCount; i++)
                {
                    double curr = grid[i][p];
                    double prev = grid[i - 1][p];

                    if(!Double.isNaN(curr) && !Double.isNaN(prev)
                    && abs(curr - prev) > RESEG_GC_SPARSE_MIN_JUMP && binProportion[i] < RESEG_GC_SPARSE_MAX_PROPORTION)
                    {
                        grid[i][p] = NaN;
                    }
                }
            }

            anyNaN = false;

            for(int p = 0; p < percentileCount; p++)
            {
                for(int i = 0; i < binCount; i++)
                {
                    if(Double.isNaN(grid[i][p]))
                        anyNaN = true;
                }
            }

            for(int p = 0; p < percentileCount; p++)
            {
                for(int i = 1; i < binCount; i++)
                {
                    if(Double.isNaN(grid[i][p]))
                        grid[i][p] = grid[i - 1][p];
                }
            }
        }
    }

    private static class RatioBafCount
    {
        public final double Ratio;
        public int BafCount;

        public RatioBafCount(double ratio)
        {
            Ratio = ratio;
            BafCount = 0;
        }

        public String toString() { return format("%.4f: %d", Ratio, BafCount); }
    }

    private static List<RatioBafCount> buildRatioBafCounts(final List<ObservedRegion> segments)
    {
        List<RatioBafCount> ratioBafCounts = Lists.newArrayList();

        for(ObservedRegion segment : segments)
        {
            RatioBafCount ratioBafCount = ratioBafCounts.stream()
                    .filter(x -> x.Ratio == segment.observedTumorRatio()).findFirst().orElse(null);

            if(ratioBafCount == null)
            {
                ratioBafCount = new RatioBafCount(segment.observedTumorRatio());
                ratioBafCounts.add(ratioBafCount);
            }

            ratioBafCount.BafCount += segment.bafCount();
        }

        Collections.sort(ratioBafCounts, Comparator.comparingDouble(x -> x.Ratio));

        return ratioBafCounts;
    }

    private static Double weightRatioAcrossPercentile(final List<RatioBafCount> ratioBafCounts, double percentile)
    {
        double totalWeight = ratioBafCounts.stream().mapToDouble(x -> x.BafCount).sum();

        if(totalWeight <= 0)
            return null;

        double cumulative = 0;
        Double lower = null;
        Double upper = null;

        for(RatioBafCount ratioBafCount : ratioBafCounts)
        {
            cumulative += ratioBafCount.BafCount;
            double cumulativeProportion = cumulative / totalWeight;

            if(lower == null && cumulativeProportion >= percentile)
                lower = ratioBafCount.Ratio;

            if(upper == null && cumulativeProportion > percentile)
                upper = ratioBafCount.Ratio;

            if(lower != null && upper != null)
                break;
        }

        if(lower == null || upper == null)
            return null;

        return (lower + upper) / 2;
    }

    private static GcBin findBin(double gcContent, final List<GcBin> allBins)
    {
        return allBins.stream().filter(x -> x.contains(gcContent)).findFirst().orElse(null);
    }
}
