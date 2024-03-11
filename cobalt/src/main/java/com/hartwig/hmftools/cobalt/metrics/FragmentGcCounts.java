package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.Math.round;
import static java.lang.String.format;

public class FragmentGcCounts implements Comparable<FragmentGcCounts>
{
    public final int FragmentLength;
    public final double GcPercent;

    public int Count;
    public int DuplicateReadCount;

    public FragmentGcCounts(final int fragmentLength, final double gcPercent)
    {
        FragmentLength = fragmentLength;
        GcPercent = gcPercent;
        Count = 0;
        DuplicateReadCount = 0;
    }

    @Override
    public int compareTo(final FragmentGcCounts other)
    {
        if(FragmentLength != other.FragmentLength)
            return FragmentLength < other.FragmentLength ? -1 : 1;

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

    public static String formKey(final int fragmentLength, final double gcPercent)
    {
        return format("%d_%.2f", fragmentLength, gcPercent);
    }
}
