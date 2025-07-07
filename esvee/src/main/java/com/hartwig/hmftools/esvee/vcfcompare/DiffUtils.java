package com.hartwig.hmftools.esvee.vcfcompare;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.List;

public final class DiffUtils
{
    private static final int DEFAULT_MAX_DIFF = 20;
    private static final double DEFAULT_MAX_DIFF_PERC = 0.2;

    private static boolean hasDiffWithinTolerance(double value1, double value2)
    {
        if(value1 <= 0 && value2 <= 0)
            return true;

        double diff = abs(value1 - value2);
        double diffPerc = diff / max(value1, value2);
        return diff <= DEFAULT_MAX_DIFF && diffPerc <= DEFAULT_MAX_DIFF_PERC;
    }

    public static void checkValueDifference(final List<String> diffs, final String type, final int value1, final int value2)
    {
        if(hasDiffWithinTolerance(value1, value2))
            return;

        diffs.add(format("%s(%d/%d)", type, value1, value2));
    }

    public static void checkValueDifference(final List<String> diffs, final String type, final double value1, final double value2)
    {
        if(hasDiffWithinTolerance(value1, value2))
            return;

        diffs.add(format("%s(%.1f/%.1f)", type, value1, value2));
    }

    public static void checkValueDifference(final List<String> diffs, final String type, final boolean value1, final boolean value2)
    {
        if(value1 == value2)
            return;

        diffs.add(format("%s(%s/%s)", type, value1, value2));
    }

    public static void checkValueDifference(final List<String> diffs, final String type, final String value1, final String value2)
    {
        if(value1.equals(value2))
            return;

        diffs.add(format("%s(%s/%s)", type, value1, value2));
    }
}
