package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

public class PartitionStats
{
    public long TotalReads;
    public int JunctionCount;
    public int JunctionFragmentCount;
    public int InitialSupportingFragmentCount;
    public int SupportingFragmentCount;
    public int LocalCompleteGroups;
    public int LocalIncompleteGroups;
    public int SpanningGroups;

    public final int[] ReadFilterCounts;

    public PartitionStats()
    {
        TotalReads = 0;
        JunctionCount = 0;
        JunctionFragmentCount = 0;
        InitialSupportingFragmentCount = 0;
        SupportingFragmentCount = 0;
        LocalCompleteGroups = 0;
        LocalIncompleteGroups = 0;
        SpanningGroups = 0;

        ReadFilterCounts = new int[ReadFilterType.values().length];
    }

    public void add(final PartitionStats other)
    {
        TotalReads += other.TotalReads;
        JunctionCount += other.JunctionCount;
        JunctionFragmentCount += other.JunctionFragmentCount;
        InitialSupportingFragmentCount += other.InitialSupportingFragmentCount;
        SupportingFragmentCount += other.SupportingFragmentCount;
        LocalCompleteGroups += other.LocalCompleteGroups;
        LocalIncompleteGroups += other.LocalIncompleteGroups;
        SpanningGroups += other.SpanningGroups;

        for(int i = 0; i < ReadFilterCounts.length; ++i)
        {
            ReadFilterCounts[i] += other.ReadFilterCounts[i];
        }
    }

    public String toString()
    {
        return format("reads(%s) junc(%d) juncFrags(%d) supportFrags(init=%d final=%d) groups(comp=%d incomp=%d span=%d)",
                TotalReads, JunctionCount, JunctionFragmentCount, InitialSupportingFragmentCount, SupportingFragmentCount,
                LocalCompleteGroups, LocalIncompleteGroups, SpanningGroups);
    }

}
