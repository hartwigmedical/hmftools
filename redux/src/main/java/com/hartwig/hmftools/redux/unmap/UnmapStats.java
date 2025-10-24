package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import java.util.concurrent.atomic.AtomicInteger;

public class UnmapStats
{
    // types of reads unmapped
    public AtomicInteger ReadCount;
    public AtomicInteger MateCount;
    public AtomicInteger SuppAlignmentCount; // supplementary alignment (ie attribute) is unmapped
    public AtomicInteger SupplementaryCount; // supplementary read is unmapped
    public AtomicInteger SecondaryCount; // secondary read is unmapped
    public AtomicInteger FullyUnmappedCount; // ie both read and mate
    public AtomicInteger ExistingUnmapped; // read and mate were already unmapped in the BAM

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
        SecondaryCount = new AtomicInteger();
        FullyUnmappedCount = new AtomicInteger();
        HighDepthCount = new AtomicInteger();
        LongSoftClipCount = new AtomicInteger();
        ChimericCount = new AtomicInteger();
        ExistingUnmapped = new AtomicInteger();
    }

    public void clearAll()
    {
        ReadCount.set(0);
        MateCount.set(0);
        SuppAlignmentCount.set(0);
        SupplementaryCount.set(0);
        SecondaryCount.set(0);
        FullyUnmappedCount.set(0);
        HighDepthCount.set(0);
        LongSoftClipCount.set(0);
        ChimericCount.set(0);
        ExistingUnmapped.set(0);
    }

    public String toString()
    {
        return format("reads(%d) mates(%d) both(%d existing=%d) supps(read=%d align=%d) seconds(%d) reasons(depth=%d softClip=%d chimeric=%d)",
                ReadCount.get(), MateCount.get(), FullyUnmappedCount.get(), ExistingUnmapped.get(), SupplementaryCount.get(), SuppAlignmentCount.get(),
                SecondaryCount.get(), HighDepthCount.get(), LongSoftClipCount.get(), ChimericCount.get());
    }

    public String unpairedStats()
    {
        return format("reads(%d) existingUnmapped(%d) supps(read=%d align=%d) seconds(%d) reasons(depth=%d softClip=%d)",
                FullyUnmappedCount.get(), ExistingUnmapped.get(), SupplementaryCount.get(), SuppAlignmentCount.get(),
                SecondaryCount.get(), HighDepthCount.get(), LongSoftClipCount.get());
    }
}
