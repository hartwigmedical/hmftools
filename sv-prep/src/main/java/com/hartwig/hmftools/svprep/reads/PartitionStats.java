package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

public class PartitionStats
{
    public long TotalReads;
    public int JunctionCount;
    public int JunctionFragmentCount;
    public int InitialSupportingFragmentCount;
    public int SupportingFragmentCount;

    public final int[] ReadFilterCounts;

    public PartitionStats()
    {
        TotalReads = 0;
        JunctionCount = 0;
        JunctionFragmentCount = 0;
        InitialSupportingFragmentCount = 0;
        SupportingFragmentCount = 0;
        ReadFilterCounts = new int[ReadFilterType.values().length];
    }

    public void add(final PartitionStats other)
    {
        TotalReads += other.TotalReads;
        JunctionCount += other.JunctionCount;
        JunctionFragmentCount += other.JunctionFragmentCount;
        InitialSupportingFragmentCount += other.InitialSupportingFragmentCount;
        SupportingFragmentCount += other.SupportingFragmentCount;
    }

    public String toString()
    {
        return format("reads(%s) junc(%d) juncFrags(%d) supportFrags(init=%d final=%d)",
                TotalReads, JunctionCount, JunctionFragmentCount, InitialSupportingFragmentCount, SupportingFragmentCount);
    }

}
