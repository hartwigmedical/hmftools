package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.Math.round;
import static java.lang.String.format;

public class FragmentGcCounts implements Comparable<FragmentGcCounts>
{
    public final int FragmentLength;
    public final int DuplicateCount;
    public final double GcPercent;

    public int Count;

    public FragmentGcCounts(final int fragmentLength, final int duplicateCount, final double gcPercent)
    {
        FragmentLength = fragmentLength;
        DuplicateCount = duplicateCount;
        GcPercent = gcPercent;
        Count = 0;
    }

    @Override
    public int compareTo(final FragmentGcCounts other)
    {
        if(FragmentLength != other.FragmentLength)
            return FragmentLength < other.FragmentLength ? -1 : 1;

        if(DuplicateCount != other.DuplicateCount)
            return DuplicateCount < other.DuplicateCount ? -1 : 1;

        if(GcPercent != other.GcPercent)
            return GcPercent < other.GcPercent ? -1 : 1;

        return 0;
    }

    public static double roundGcPercent(double gcPercent, double gcUnits)
    {
        return gcUnits * round(gcPercent/ gcUnits);
    }

    public static int roundFragmentLength(int fragmentLength, int lengthUnits)
    {
        return (int)(lengthUnits * round(fragmentLength /(double)lengthUnits));
    }

    private static final double DUP_ROUNDING = 10.0;
    private static final int DUP_ROUNDING_THRESHOLD = 20;

    public static int roundDuplicateCount(final int duplicateCount)
    {
        if(duplicateCount < DUP_ROUNDING_THRESHOLD)
            return duplicateCount;

        return (int)(DUP_ROUNDING * round(duplicateCount / DUP_ROUNDING));
    }

    public static String formKey(final int fragmentLength, final int duplicateCount, final double gcPercent)
    {
        return format("%d_%d_%.2f", fragmentLength, duplicateCount, gcPercent);
    }
}
