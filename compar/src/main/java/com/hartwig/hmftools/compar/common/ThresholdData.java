package com.hartwig.hmftools.compar.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

public class ThresholdData
{
    public final ThresholdType Type;
    public final double AbsoluteDiff;
    public final double PercentDiff;

    public ThresholdData(final ThresholdType type, final double absoluteDiff, final double percentDiff)
    {
        Type = type;
        AbsoluteDiff = absoluteDiff;
        PercentDiff = percentDiff;
    }

    public boolean hasDiff(double value1, double value2)
    {
        if(value1 == 0 && value2 == 0)
            return false;

        double absDiff = abs(value1 - value2);

        boolean hasAbsDiff = absDiff > AbsoluteDiff;
        boolean hasRelDiff = absDiff / max(abs(value1), abs(value2)) > PercentDiff;

        if(Type == ThresholdType.ABSOLUTE_AND_PERCENT)
            return hasAbsDiff && hasRelDiff;

        if(Type == ThresholdType.ABSOLUTE)
            return hasAbsDiff;

        return hasRelDiff;
    }

    public String toString() { return format("type(%s) abs(%f) perc(%.3f)", Type, AbsoluteDiff, PercentDiff); }
}
