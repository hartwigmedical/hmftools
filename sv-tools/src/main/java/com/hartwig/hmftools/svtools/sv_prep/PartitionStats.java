package com.hartwig.hmftools.svtools.sv_prep;

import static java.lang.String.format;

public class PartitionStats
{
    public long TotalReads;
    public int Buckets;
    public int EmptyBuckets;
    public int JunctionCount;
    public int JunctionFragmentCount;
    public int InitialSupportingReadCount;
    public int SupportingReadCount;

    public final int[] ReadFilterCounts;

    public PartitionStats()
    {
        TotalReads = 0;
        Buckets = 0;
        EmptyBuckets = 0;
        JunctionCount = 0;
        JunctionFragmentCount = 0;
        InitialSupportingReadCount = 0;
        SupportingReadCount = 0;
        ReadFilterCounts = new int[ReadFilterType.values().length];
    }

    public void add(final PartitionStats other)
    {
        TotalReads += other.TotalReads;
        Buckets += other.Buckets;
        EmptyBuckets += other.EmptyBuckets;
        JunctionCount += other.JunctionCount;
        JunctionFragmentCount += other.JunctionFragmentCount;
        InitialSupportingReadCount += other.InitialSupportingReadCount;
        SupportingReadCount += other.SupportingReadCount;
    }

    public String toString()
    {
        return format("reads(%s) buckets(%d empty=%d) junc(%d) juncFrags(%d) support(init=%d final=%d)",
                TotalReads, Buckets, EmptyBuckets, JunctionCount, JunctionFragmentCount, InitialSupportingReadCount, SupportingReadCount);
    }

}
