package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.Arrays;
import java.util.List;

// Distribution summary over a set of integer values.
public record SummaryStats(double mean, double min, double p5, double p50, double p95, double max)
{
    static SummaryStats from(int[] values)
    {
        int[] sorted = values.clone();
        Arrays.sort(sorted);
        return fromSorted(sorted);
    }

    static SummaryStats from(List<Integer> values)
    {
        return fromSorted(values.stream().mapToInt(Integer::intValue).sorted().toArray());
    }

    private static SummaryStats fromSorted(int[] sorted)
    {
        if(sorted.length == 0)
        {
            throw new IllegalArgumentException("cannot summarise an empty set of values");
        }

        long sum = 0;
        for(int value : sorted)
        {
            sum += value;
        }
        double mean = (double) sum / sorted.length;

        return new SummaryStats(
                mean, sorted[0], percentile(sorted, 5), percentile(sorted, 50), percentile(sorted, 95), sorted[sorted.length - 1]);
    }

    // Nearest-rank percentile over an already-sorted array.
    private static double percentile(int[] sorted, double percent)
    {
        int index = (int) ceil(percent / 100.0 * sorted.length) - 1;
        return sorted[max(0, min(sorted.length - 1, index))];
    }
}
