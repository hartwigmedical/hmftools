package com.hartwig.hmftools.markdups.consensus;

import static java.lang.String.format;

public class ConsensusStatistics
{
    private int mDualStrandMismatchReadGroupCount;
    private int mDualStrandMismatchReadCount;

    public ConsensusStatistics()
    {
        mDualStrandMismatchReadGroupCount = 0;
        mDualStrandMismatchReadCount = 0;
    }

    public void registerDualStrandMismatchReadGroup(int readCount)
    {
        ++mDualStrandMismatchReadGroupCount;
        mDualStrandMismatchReadCount += readCount;
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