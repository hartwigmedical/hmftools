package com.hartwig.hmftools.common.utils;

import static java.lang.Math.pow;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public final class Doubles
{
    private static final double EPSILON = 1e-10;

    public static boolean equal(double first, double second) {
        return Math.abs(first - second) < EPSILON;
    }

    public static boolean isZero(double value) {
        return equal(value, 0);
    }

    public static boolean lessThan(double value, double reference) {
        return value - reference < -EPSILON;
    }

    public static boolean lessOrEqual(double value, double reference) {
        return value - reference < EPSILON;
    }

    public static boolean greaterThan(double value, double reference) {
        return value - reference > EPSILON;
    }

    public static boolean greaterOrEqual(double value, double reference) {
        return value - reference > -EPSILON;
    }

    public static double replaceNaNWithZero(double value) { return Double.isNaN(value) ? 0d : value; }

    public static boolean positive(double value) {
        return greaterThan(value, 0);
    }

    public static boolean positiveOrZero(double value) {
        return greaterOrEqual(value, 0);
    }

    public static double absDistanceFromInteger(double value) {
        return Math.abs((value - Math.round(value)));
    }

    public static Comparator<Double> comparator() {
        return (o1, o2) -> Doubles.equal(o1, o2) ? 0 :  Double.compare(o1, o2);
    }

    public static double median(final List<Double> values) {
        return median(values, x -> true);
    }

    public static double median(final List<Double> values, @NotNull final Predicate<Double> filter)
    {
        final List<Double> reads = values.stream().filter(filter).sorted().collect(Collectors.toList());
        return medianOfSortedDoubles(reads);
    }

    private static double medianOfSortedDoubles(final List<Double> values)
    {
        int count = values.size();

        if (count == 0)
            return 0;

        return count % 2 == 0 ? (values.get(count / 2) + values.get(count / 2 - 1)) / 2 : values.get(count / 2);
    }

    public static int medianInteger(final List<Integer> values)
    {
        Collections.sort(values);

        int count = values.size();

        if(count == 0)
            return 0;

        return count % 2 == 0 ? (values.get(count / 2) + values.get(count / 2 - 1)) / 2 : values.get(count / 2);
    }

    public static double round(final double value, final int decimalPlaces)
    {
        double logValue = Math.round(pow(10, decimalPlaces));
        return Math.round(value * logValue) / logValue;
    }

    public static double interpolatedMedian(Collection<Double> input)
    {
        if (input.isEmpty())
        {
            return 0;
        }
        List<Double> values = input.stream().sorted().collect(Collectors.toList());
        int count = values.size();

        double median = medianOfSortedDoubles(values);

        // now count how many are below, how many are above
        int countBelowMed = 0;
        int countAboveMed = 0;

        int cumulativeCount = 0;

        for(double val : values)
        {
            if(val != (int)val)
            {
                throw new IllegalArgumentException("inputs to interpolatedMedian must be integer values");
            }

            if(cumulativeCount == count / 2.0 && val != median)
            {
                // same corner case as median, if it falls between two numbers, say 2 and 3, then
                // we return 2.5
                return median;
            }

            if(Doubles.lessThan(val, median))
            {
                countBelowMed++;
            }
            else if(Doubles.greaterThan(val, median))
            {
                countAboveMed++;
            }

            cumulativeCount++;
        }

        double l = ((double)countBelowMed) / count;
        double r = ((double)countAboveMed) / count;

        return median - 0.5 + (0.5 - l) / (1 - l - r);
    }
}
