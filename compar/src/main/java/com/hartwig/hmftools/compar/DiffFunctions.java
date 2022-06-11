package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public final class DiffFunctions
{
    public static boolean checkDiff(final List<String> diffs, final String name, int value1, int value2, final DiffThresholds thresholds)
    {
        if(thresholds.hasDifference(name, value1, value2))
        {
            diffs.add(String.format("%s(%d/%d)", name, value1, value2));
            return true;
        }

        return false;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, double value1, double value2, final DiffThresholds thresholds)
    {
        if(thresholds.hasDifference(name, value1, value2))
        {
            diffs.add(String.format("%s(%.3f/%.3f)", name, value1, value2));
            return true;
        }

        return false;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, int value1, int value2)
    {
        if(value1 != value2)
        {
            diffs.add(String.format("%s(%d/%d)", name, value1, value2));
            return true;
        }

        return false;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, boolean value1, boolean value2)
    {
        if(value1 == value2)
            return false;

        diffs.add(String.format("%s(%s/%s)", name, value1, value2));
        return true;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, String value1, final String value2)
    {
        if(value1.equals(value2))
            return false;

        diffs.add(String.format("%s(%s/%s)", name, value1, value2));
        return true;
    }

    public static void checkFilterDiffs(Set<String> refFilters, final Set<String> otherFilters, final List<String> diffs)
    {
        Set<String> origFilterDiffs = refFilters.stream().filter(x -> !otherFilters.contains(x)).collect(Collectors.toSet());
        Set<String> newFilterDiffs = otherFilters.stream().filter(x -> !refFilters.contains(x)).collect(Collectors.toSet());

        if(!newFilterDiffs.isEmpty() || !origFilterDiffs.isEmpty())
        {
            diffs.add(String.format("%s(%s/%s)", "filter", filtersStr(origFilterDiffs), filtersStr(newFilterDiffs)));
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
