package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

public class FragmentLengthBounds
{
    public final int UpperBound;
    public final int LowerBound;
    public final int Median;
    public final double StdDeviation;

    public static final FragmentLengthBounds INVALID = new FragmentLengthBounds(0, 0, 0, 0);

    public FragmentLengthBounds(final int lowerBound, final int upperBound, final int median, final double stdDeviation)
    {
        UpperBound = upperBound;
        LowerBound = lowerBound;
        Median = median;
        StdDeviation = stdDeviation;
    }

    public boolean isValid() { return LowerBound > 0 && UpperBound > LowerBound; }

    public String toString() { return format("bounds(%d - %d) median(%d stdDev=%.1f)",
            LowerBound, UpperBound, Median, StdDeviation); }
}
