package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ResolvedFragmentState
{
    public final FragmentStatus Status;
    public final String Coordinates;

    public boolean MateReceived;
    public int ExpectedSupplementaries;
    public int ProcessedSupplementaries;

    public ResolvedFragmentState(
            final FragmentStatus status, final String coordinates,
            final int expectedSupplementaries, final int processedSupplementaries, final boolean mateReceived)
    {
        Status = status;
        Coordinates = coordinates;
        MateReceived = mateReceived;
        ExpectedSupplementaries = expectedSupplementaries;
        ProcessedSupplementaries = processedSupplementaries;
    }

    public boolean allReceived()
    {
        return MateReceived && ExpectedSupplementaries == ProcessedSupplementaries;
    }

    public static ResolvedFragmentState fragmentState(final Fragment fragment)
    {
        int nonSuppCount = 0;
        int expectedSuppCount = 0;
        int processedSuppCount = 0;

        for(SAMRecord read : fragment.reads())
        {
            if(read.getSupplementaryAlignmentFlag())
            {
                ++processedSuppCount;
            }
            else
            {
                ++nonSuppCount;

                expectedSuppCount += SupplementaryReadData.alignmentCount(read);
            }
        }

        return new ResolvedFragmentState(
                fragment.status(), fragment.coordinates().Key,
                expectedSuppCount, processedSuppCount, nonSuppCount > 1 || fragment.unpaired());
    }

    public void update(final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag())
        {
            ++ProcessedSupplementaries;
        }
        else
        {
            MateReceived = true;
            ExpectedSupplementaries += SupplementaryReadData.alignmentCount(read);
        }
    }

    public String toString() { return format("status(%s) coords(%s) mate(%s) supps(%d/%d)",
            Status, Coordinates, MateReceived ? "received" : "pending", ProcessedSupplementaries, ExpectedSupplementaries); }
}
