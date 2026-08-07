package com.hartwig.hmftools.purple.reseg;

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
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class GcNormaliser
{
    // per-sample bucketed GC-content normalisation of observedTumorRatio, restricted to diploid
    // segments with a positive bafCount (and, if a primary ratio peak was identified, within its bounds)
    private GcNormaliser() {}

    private static final int MAX_SUPPRESSION_PASSES = 20;

    public static List<ObservedRegion> normalise(final List<ObservedRegion> allSegments, final Optional<RatioPeakResult> peakResult)
    {
        List<GcBin> allBins = GcBin.allBins();

        // the set of bins actually touched by any segment in the full sample, used to resolve the
        // out-of-range LOW/HIGH clamp bins to the observed extremes rather than the theoretical ones
        Set<GcBin> observedBins = Sets.newHashSet();
        for(ObservedRegion segment : allSegments)
        {
            findBin(segment.gcContent(), allBins).ifPresent(observedBins::add);
        }

        GcBin lowClampBin = observedBins.isEmpty()
                ? allBins.get(0) : Collections.min(observedBins, Comparator.comparingDouble(b -> b.Lower));

        GcBin highClampBin = observedBins.isEmpty()
                ? allBins.get(allBins.size() - 1) : Collections.max(observedBins, Comparator.comparingDouble(b -> b.Upper));

        List<ObservedRegion> population = allSegments.stream()
                .filter(x -> x.germlineStatus() == GermlineStatus.DIPLOID && x.bafCount() > 0)
                .filter(x -> peakResult.map(r -> x.observedTumorRatio() >= r.LeftBound && x.observedTumorRatio() <= r.RightBound).orElse(true))
                .collect(Collectors.toList());

        Map<GcBin, Double> gcAdjustments = calculateGcAdjustments(population, allBins);

        List<ObservedRegion> result = new ArrayList<>(allSegments.size());

        for(ObservedRegion segment : allSegments)
        {
            GcBin bin = findBin(segment.gcContent(), allBins)
                    .orElseGet(() -> segment.gcContent() <= lowClampBin.Lower ? lowClampBin : highClampBin);

            double adjustment = gcAdjustments.getOrDefault(bin, 1.0);

            if(Double.isNaN(adjustment) || adjustment <= 0)
                adjustment = 1.0;

            double newRatio = round(segment.observedTumorRatio() / adjustment, 4);

            ObservedRegion newRegion = ObservedRegion.fromOther(segment);
            newRegion.setObservedTumorRatio(newRatio);
            result.add(newRegion);
        }

        return result;
    }

    static Map<GcBin, Double> calculateGcAdjustments(final List<ObservedRegion> population, final List<GcBin> allBins)
    {
        double totalWeight = population.stream().mapToDouble(s -> s.bafCount()).sum();

        double[] binProportion = new double[allBins.size()];

        for(int i = 0; i < allBins.size(); i++)
        {
            GcBin bin = allBins.get(i);
            double binWeight = population.stream().filter(x -> bin.contains(x.gcContent())).mapToDouble(s -> s.bafCount()).sum();
            binProportion[i] = totalWeight > 0 ? binWeight / totalWeight : 0;
        }

        double[][] grid = new double[allBins.size()][RESEG_GC_PERCENTILES.length];

        for(int p = 0; p < RESEG_GC_PERCENTILES.length; p++)
        {
            double percentile = RESEG_GC_PERCENTILES[p];

            OptionalDouble allValue = weightedPercentileAcross(population, s -> true, percentile);

            for(int i = 0; i < allBins.size(); i++)
            {
                GcBin bin = allBins.get(i);
                OptionalDouble binValue = weightedPercentileAcross(population, x -> bin.contains(x.gcContent()), percentile);

                grid[i][p] = (binValue.isPresent() && allValue.isPresent() && allValue.getAsDouble() != 0)
                        ? binValue.getAsDouble() / allValue.getAsDouble()
                        : Double.NaN;
            }
        }

        suppressSparseValues(grid, binProportion);

        Map<GcBin, Double> adjustments = new LinkedHashMap<>();

        for(int i = 0; i < allBins.size(); i++)
        {
            List<Double> valid = new ArrayList<>();

            for(double v : grid[i])
            {
                if(!Double.isNaN(v))
                    valid.add(v);
            }

            adjustments.put(allBins.get(i), valid.isEmpty() ? Double.NaN : median(valid));
        }

        return adjustments;
    }

    // marks (percentile,bin) entries that jump too much from the preceding bin and are backed by too
    // little data as unreliable, then forward-fills them from the nearest preceding valid bin, repeating
    // until a full pass finds nothing left to suppress
    private static void suppressSparseValues(final double[][] grid, final double[] binProportion)
    {
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
                    && Math.abs(curr - prev) > RESEG_GC_SPARSE_MIN_JUMP && binProportion[i] < RESEG_GC_SPARSE_MAX_PROPORTION)
                    {
                        grid[i][p] = Double.NaN;
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

    private static OptionalDouble weightedPercentileAcross(
            final List<ObservedRegion> segments, final Predicate<ObservedRegion> filter, double percentile)
    {
        Map<Double, Double> weightByRatio = new TreeMap<>();

        for(ObservedRegion segment : segments)
        {
            if(!filter.test(segment))
                continue;

            weightByRatio.merge(segment.observedTumorRatio(), (double)segment.bafCount(), Double::sum);
        }

        List<WeightedPercentile.RatioWeightPair> pairs = weightByRatio.entrySet().stream()
                .map(e -> new WeightedPercentile.RatioWeightPair(e.getKey(), e.getValue()))
                .collect(Collectors.toList());

        return WeightedPercentile.compute(pairs, percentile);
    }

    private static Optional<GcBin> findBin(double gcContent, final List<GcBin> allBins)
    {
        for(GcBin bin : allBins)
        {
            if(bin.contains(gcContent))
                return Optional.of(bin);
        }

        return Optional.empty();
    }
}
