package com.hartwig.hmftools.markdups.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.List;

import htsjdk.samtools.SAMRecord;

public class ConsensusStatistics
{
    private int mDualStrandMismatchReadGroupCount;
    private int mDualStrandMismatchReadCount;

    public ConsensusStatistics()
    {
        mDualStrandMismatchReadGroupCount = 0;
        mDualStrandMismatchReadCount = 0;
    }

    public void registerDualStrandMismatchReadGroup(final List<SAMRecord> reads)
    {
        MD_LOGGER.trace("Mismatched basis found in dual stranded read group readCount({})", reads.size());
        for(SAMRecord read : reads)
        {
            MD_LOGGER.trace("Read in dual stranded read group: {}", read);
        }

        ++mDualStrandMismatchReadGroupCount;
        mDualStrandMismatchReadCount += reads.size();
    }

    public void merge(final ConsensusStatistics other)
    {
        mDualStrandMismatchReadGroupCount += other.mDualStrandMismatchReadGroupCount;
        mDualStrandMismatchReadCount += other.mDualStrandMismatchReadCount;
    }

    @Override
    public String toString()
    {
        return format("dualStrandMismatchReadGroupCount(%d) dualStrandMismatchReadCount(%d)", mDualStrandMismatchReadGroupCount, mDualStrandMismatchReadCount);
    }
}