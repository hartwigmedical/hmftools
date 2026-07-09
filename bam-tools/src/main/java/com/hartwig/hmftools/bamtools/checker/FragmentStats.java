package com.hartwig.hmftools.bamtools.checker;

import static java.lang.Math.max;
import static java.lang.String.format;

public class FragmentStats
{
    public long TotalFragments;
    public long InterPartitionFragments;
    public long FragmentsWithSupplementaries;
    public long MaxFragmentCount;
    public long MateCigarFixed;
    public long SuppsDropped;
    public long PrimariesUnmapped;
    public long NonRefContigDropped;

    public FragmentStats()
    {
        TotalFragments = 0;
        InterPartitionFragments = 0;
        FragmentsWithSupplementaries = 0;
        MaxFragmentCount = 0;
        MateCigarFixed = 0;
        SuppsDropped = 0;
        PrimariesUnmapped = 0;
        NonRefContigDropped = 0;
    }

    public void reset()
    {
        TotalFragments = 0;
        InterPartitionFragments = 0;
        FragmentsWithSupplementaries = 0;
        MaxFragmentCount = 0;
        MateCigarFixed = 0;
        SuppsDropped = 0;
        PrimariesUnmapped = 0;
        NonRefContigDropped = 0;
    }

    public void merge(final FragmentStats other)
    {
        TotalFragments += other.TotalFragments;
        InterPartitionFragments += other.InterPartitionFragments;
        FragmentsWithSupplementaries += other.FragmentsWithSupplementaries;
        MaxFragmentCount = max(MaxFragmentCount, other.MaxFragmentCount);
        MateCigarFixed += other.MateCigarFixed;
        PrimariesUnmapped += other.PrimariesUnmapped;
        SuppsDropped += other.SuppsDropped;
        NonRefContigDropped += other.NonRefContigDropped;
    }

    public String toString()
    {
        return format("fragments(%d) interPartition(%d) withSupp(%d) maxFrags(%d) mateCigarFixed(%d) primaryUnmapped(%d) suppDropped(%d) nonRefContigDropped(%d)",
                TotalFragments, InterPartitionFragments, FragmentsWithSupplementaries, MaxFragmentCount, MateCigarFixed,
                PrimariesUnmapped, SuppsDropped, NonRefContigDropped);
    }
}
