package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;

import htsjdk.samtools.SAMRecord;

class ResolvedFragmentState
{
    public final FragmentStatus Status;

    public boolean MateReceived;
    public int ExpectedSupplementaries;
    public int ProcessedSupplementaries;

    public ResolvedFragmentState(final FragmentStatus status, final int expectedSupplementaries, final int processedSupplementaries, final boolean mateReceived)
    {
        Status = status;
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

                if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
                {
                    ++expectedSuppCount;
                }
            }
        }

        return new ResolvedFragmentState(
                fragment.status(), expectedSuppCount, processedSuppCount, nonSuppCount > 1 || fragment.unpaired());
    }

    public void update(final List<SAMRecord> reads)
    {
        for(SAMRecord read : reads)
        {
            if(read.getSupplementaryAlignmentFlag())
            {
                ++ProcessedSupplementaries;
            }
            else
            {
                MateReceived = true;

                if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
                {
                    ++ExpectedSupplementaries;
                }
            }
        }
    }

    public boolean isValid()
    {
        // can receive a supplementary for a mate before the mate is received
        return ProcessedSupplementaries <= 2 && ExpectedSupplementaries <= 2;
        // return ProcessedSupplementaries <= ExpectedSupplementaries && ExpectedSupplementaries <= 2;
    }

    public String toString() { return format("status(%s) mate(%s) supps(%d/%d)",
            Status, MateReceived ? "received" : "pending", ProcessedSupplementaries, ExpectedSupplementaries); }
}
