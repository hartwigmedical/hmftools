package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

public class ReadCounts
{
    public long Total;
    public long Duplicates;
    public long DualStrand;

    public ReadCounts()
    {
        Total = 0;
        Duplicates = 0;
        DualStrand = 0;
    }

    public void merge(final ReadCounts other)
    {
        Total += other.Total;
        Duplicates += other.Duplicates;
        DualStrand += other.DualStrand;
    }

    public String toString()
    {
        return format("total=%d dup=%d dual=%d", Total, Duplicates, DualStrand);
    }
}
