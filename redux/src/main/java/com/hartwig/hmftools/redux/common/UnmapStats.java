package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.concurrent.atomic.AtomicInteger;

public class UnmapStats
{
    // types of reads unmapped
    public AtomicInteger ReadCount;
    public AtomicInteger MateCount;
    public AtomicInteger SuppAlignmentCount; // supplementary alignment (ie attribute) is unmapped
    public AtomicInteger SupplementaryCount; // supplementary read is unmapped
    public AtomicInteger UnmappedCount; // ie both read and mate

    // reasons for unmapping a read
    public AtomicInteger HighDepthCount;
    public AtomicInteger LongSoftClipCount;
    public AtomicInteger ChimericCount;

    public UnmapStats()
    {
        ReadCount = new AtomicInteger();
        MateCount = new AtomicInteger();
        SuppAlignmentCount = new AtomicInteger();
        SupplementaryCount = new AtomicInteger();
        UnmappedCount = new AtomicInteger();
        HighDepthCount = new AtomicInteger();
        LongSoftClipCount = new AtomicInteger();
        ChimericCount = new AtomicInteger();
    }

    public String toString()
    {
        return format("reads(%d) mates(%d) both(%d) supplementary(read=%d align=%d) reasons(depth=%d softClip=%d chimeric=%d)",
            ReadCount.get(), MateCount.get(), UnmappedCount.get(), SupplementaryCount.get(), SuppAlignmentCount.get(),
                HighDepthCount.get(), LongSoftClipCount.get(), ChimericCount.get());
    }
}
