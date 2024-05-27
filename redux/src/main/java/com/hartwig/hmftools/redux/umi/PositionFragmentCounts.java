package com.hartwig.hmftools.redux.umi;

import static java.lang.String.format;

public class PositionFragmentCounts
{
    // counts are keyed by the first 2 elements:
    public final int UniqueCoordCount;
    public final int UniqueFragmentCount;

    public int Frequency;
    public int MaxCoordUmiCount;
    public int MaxUmiReadsCount; // max (primary) reads collapsed into a UMI group
    public String UmiGroupDetails;

    public PositionFragmentCounts(final int uniqueCoordCount, final int uniqueFragmentCount)
    {
        UniqueCoordCount = uniqueCoordCount;
        UniqueFragmentCount = uniqueFragmentCount;
        Frequency = 0;
        MaxCoordUmiCount = 0;
        MaxUmiReadsCount = 0;
        UmiGroupDetails = "";
    }

    public String toString()
    {
        return format("pos(%d) unique(%d) count(%d) max(umi=%d reads=%d)",
                UniqueCoordCount, UniqueFragmentCount, Frequency, MaxCoordUmiCount, MaxUmiReadsCount);
    }
}
