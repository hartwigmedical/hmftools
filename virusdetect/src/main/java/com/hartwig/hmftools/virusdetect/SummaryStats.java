package com.hartwig.hmftools.virusdetect;

import java.lang.reflect.RecordComponent;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

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
    private static final List<RecordComponent> FIELDS = List.of(SummaryStats.class.getRecordComponents());

    public static final List<String> FIELD_NAMES = FIELDS.stream().map(RecordComponent::getName).toList();

    // Values by field name, in the same order.
    public Map<String, Double> fieldValues()
    {
        // TODO: stream operation. the exception handling is completely unnecessary.
        Map<String, Double> values = new LinkedHashMap<>();
        for(RecordComponent field : FIELDS)
        {
            try
            {
                values.put(field.getName(), (Double) field.getAccessor().invoke(this));
            }
            catch(ReflectiveOperationException e)
            {
                throw new IllegalStateException("Cannot read summary stats field: " + field.getName(), e);
            }
        }
        return values;
    }

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
