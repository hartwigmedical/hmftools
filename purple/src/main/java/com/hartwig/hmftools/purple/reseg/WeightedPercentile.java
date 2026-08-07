package com.hartwig.hmftools.purple.reseg;

import java.util.List;
import java.util.OptionalDouble;

public final class WeightedPercentile
{
    // observedTumorRatio at a given weighted percentile across a set of (ratio, bafCount) pairs
    // the average of the ratio values straddling the percentile crossing of cumulative weight.
    private WeightedPercentile() {}

    public record RatioWeightPair(double ratio, double weight) {}

    // sortedRatioWeights must be sorted ascending by ratio, one entry per distinct ratio value
    public static OptionalDouble compute(final List<RatioWeightPair> sortedRatioWeights, double percentile)
    {
        double totalWeight = sortedRatioWeights.stream().mapToDouble(RatioWeightPair::weight).sum();

        if(totalWeight <= 0)
            return OptionalDouble.empty();

        double cumulative = 0;
        Double lower = null;
        Double upper = null;

        for(RatioWeightPair pair : sortedRatioWeights)
        {
            cumulative += pair.weight();
            double cumulativeProportion = cumulative / totalWeight;

            if(lower == null && cumulativeProportion >= percentile)
                lower = pair.ratio();

            if(upper == null && cumulativeProportion > percentile)
                upper = pair.ratio();

            if(lower != null && upper != null)
                break;
        }

        if(lower == null || upper == null)
            return OptionalDouble.empty();

        return OptionalDouble.of((lower + upper) / 2);
    }
}
