package com.hartwig.hmftools.markdups.umi;

import static java.lang.String.format;

public class PositionFragmentCounts
{
    public final int PosGroupCount;
    public final int UniqueFragmentCount;
    public int Frequency;
    public int MaxPosUmiCount;
    public int MaxUmiReadsCount; // max reads collapsed into a UMI group
    public String UmiGroupDetails;

    public PositionFragmentCounts(final int posGroupCount, final int uniqueFragmentCount)
    {
        PosGroupCount = posGroupCount;
        UniqueFragmentCount = uniqueFragmentCount;
        Frequency = 0;
        MaxPosUmiCount = 0;
        MaxUmiReadsCount = 0;
        UmiGroupDetails = "";
    }

    public String toString()
    {
        return format("pos(%d) unique(%d) count(%d) max(umi=%d reads=%d)",
                PosGroupCount, UniqueFragmentCount, Frequency, MaxPosUmiCount, MaxUmiReadsCount);
    }
}
