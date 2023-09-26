package com.hartwig.hmftools.markdups.common;

import static java.lang.String.format;

import java.util.concurrent.atomic.AtomicInteger;

public class UnmapStats
{
    public AtomicInteger ReadCount;
    public AtomicInteger MateCount;
    public AtomicInteger SuppAlignmentCount; // supplementary alignment (ie attribute) is unmapped
    public AtomicInteger SupplementaryCount; // supplementary read is unmapped
    public AtomicInteger UnmappedCount; // ie both read and mate

    public UnmapStats()
    {
        ReadCount = new AtomicInteger();
        MateCount = new AtomicInteger();
        SuppAlignmentCount = new AtomicInteger();
        SupplementaryCount = new AtomicInteger();
        UnmappedCount = new AtomicInteger();
    }

    public String toString() { return format("reads(%d) mates(%d) both(%d) supplementary(read=%d align=%d)",
            ReadCount.get(), MateCount.get(), UnmappedCount.get(), SupplementaryCount.get(), SuppAlignmentCount.get()); }
}
