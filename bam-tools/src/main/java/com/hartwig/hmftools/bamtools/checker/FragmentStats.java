package com.hartwig.hmftools.bamtools.checker;

import static java.lang.Math.max;
import static java.lang.String.format;

public class FragmentStats
{
    public long TotalFragments;
    public long InterPartitionFragments;
    public long FragmentsWithSupplementaries;
    public long MaxFragmentCount;

    public FragmentStats()
    {
        TotalFragments = 0;
        InterPartitionFragments = 0;
        FragmentsWithSupplementaries = 0;
        MaxFragmentCount = 0;
    }

    public void reset()
    {
        TotalFragments = 0;
        InterPartitionFragments = 0;
        FragmentsWithSupplementaries = 0;
        MaxFragmentCount = 0;
    }

    public void merge(final FragmentStats other)
    {
        TotalFragments += other.TotalFragments;
        InterPartitionFragments += other.InterPartitionFragments;
        FragmentsWithSupplementaries += other.FragmentsWithSupplementaries;
        MaxFragmentCount = max(MaxFragmentCount, other.MaxFragmentCount);
    }

    public String toString()
    {
        return format("fragments(%d) interPartition(%d) withSupp(%d) maxFrags(%d)",
                TotalFragments, InterPartitionFragments, FragmentsWithSupplementaries, MaxFragmentCount);
    }
}
