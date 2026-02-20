package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.common.DiffThresholds.DEFAULT_DECIMAL_THRESHOLD;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public final class DiffFunctions
{
    public static final String FILTER_DIFF = "filter";

    public static boolean checkDiff(final List<String> diffs, final String name, long value1, long value2, final DiffThresholds thresholds)
    {
        if(!thresholds.isFieldRegistered(name))
            return checkDiff(diffs, name, value1, value2);

        if(thresholds.hasDifference(name, value1, value2))
        {
            diffs.add(format("%s(%d/%d)", name, value1, value2));
            return true;
        }

        return false;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, double value1, double value2, final DiffThresholds thresholds)
    {
        if((!thresholds.isFieldRegistered(name) && DEFAULT_DECIMAL_THRESHOLD.hasDiff(value1, value2))
        || thresholds.hasDifference(name, value1, value2))
        {
            diffs.add(format("%s(%.3f/%.3f)", name, value1, value2));
            return true;
        }

        return false;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, long value1, long value2)
    {
        if(value1 != value2)
        {
            diffs.add(format("%s(%d/%d)", name, value1, value2));
            return true;
        }

        return false;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, boolean value1, boolean value2)
    {
        if(value1 == value2)
            return false;

        diffs.add(format("%s(%s/%s)", name, value1, value2));
        return true;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, String value1, final String value2)
    {
        if(value1.equals(value2))
            return false;

        diffs.add(format("%s(%s/%s)", name, value1, value2));
        return true;
    }

    public static void checkFilterDiffs(final Set<String> refFilters, final Set<String> otherFilters, final List<String> diffs)
    {
        Set<String> refFilterDiffs = refFilters.stream().filter(x -> !otherFilters.contains(x)).collect(Collectors.toSet());
        Set<String> newFilterDiffs = otherFilters.stream().filter(x -> !refFilters.contains(x)).collect(Collectors.toSet());

        if(!newFilterDiffs.isEmpty() || !refFilterDiffs.isEmpty())
        {
            diffs.add(format("%s(%s/%s)", FILTER_DIFF, filtersStr(refFilterDiffs), filtersStr(newFilterDiffs)));
        }
    }

    public static String filtersStr(final Set<String> filters)
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        filters.forEach(x -> sj.add(x));
        return sj.toString();
    }

    public static String diffsStr(final List<String> diffs)
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        diffs.forEach(x -> sj.add(x));
        return sj.toString();
    }
}
