package com.hartwig.hmftools.viridian.common;

import static java.util.Map.entry;
import static java.util.stream.Collectors.toMap;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

// Distribution summary over a set of integer values.
public record SummaryStats(
        double mean,
        double min,
        double p5,
        double p25,
        double p50,
        double p75,
        double p95,
        double max
)
{
    // Used for writing output to reduce verbosity.
    private static final List<Map.Entry<String, ToDoubleFunction<SummaryStats>>> FIELDS = List.of(
            entry("mean", SummaryStats::mean),
            entry("min", SummaryStats::min),
            entry("p5", SummaryStats::p5),
            entry("p25", SummaryStats::p25),
            entry("p50", SummaryStats::p50),
            entry("p75", SummaryStats::p75),
            entry("p95", SummaryStats::p95),
            entry("max", SummaryStats::max));

    public static final List<String> FIELD_NAMES = FIELDS.stream().map(Map.Entry::getKey).toList();

    // Values by field name, in the same order.
    public Map<String, Double> fieldValues()
    {
        return FIELDS.stream().collect(toMap(
                Map.Entry::getKey, field -> field.getValue().applyAsDouble(this), (first, second) -> first, LinkedHashMap::new));
    }

    public static SummaryStats from(int[] values)
    {
        int[] sorted = values.clone();
        Arrays.sort(sorted);
        return fromSorted(sorted);
    }

    public static SummaryStats from(List<Integer> values)
    {
        return fromSorted(values.stream().mapToInt(Integer::intValue).sorted().toArray());
    }

    private static SummaryStats fromSorted(int[] sorted)
    {
        if(sorted.length == 0)
        {
            throw new IllegalArgumentException("Cannot summarise an empty set of values");
        }

        long sum = 0;
        for(int value : sorted)
        {
            sum += value;
        }
        double mean = (double) sum / sorted.length;

        return new SummaryStats(
                mean, sorted[0],
                percentile(sorted, 5), percentile(sorted, 25), percentile(sorted, 50), percentile(sorted, 75), percentile(sorted, 95),
                sorted[sorted.length - 1]);
    }

    // Nearest-rank percentile over a sorted array.
    private static double percentile(int[] sorted, double percent)
    {
        int index = (int) Math.ceil(percent / 100.0 * sorted.length) - 1;
        return sorted[Math.max(0, Math.min(sorted.length - 1, index))];
    }
}
