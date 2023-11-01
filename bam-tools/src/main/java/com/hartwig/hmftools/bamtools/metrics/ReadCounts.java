package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

public class ReadCounts
{
    public int TotalReads;
    public int Duplicates;
    public int DualStrand;

    public ReadCounts()
    {
        TotalReads = 0;
        Duplicates = 0;
        DualStrand = 0;
    }

    public void merge(final ReadCounts other)
    {
        TotalReads += other.TotalReads;
        Duplicates += other.Duplicates;
        DualStrand += other.DualStrand;
    }

    public String toString()
    {
        return format("total=%d dup=%d dual=%d", TotalReads, Duplicates, DualStrand);
    }
}
