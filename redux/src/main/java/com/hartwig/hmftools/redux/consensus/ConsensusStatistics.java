package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import java.util.StringJoiner;

public class ConsensusStatistics
{
    private int mDualStrandMismatchReadGroupCount;
    private int mDualStrandMismatchReadCount;
    private final int[] mOutcomeCounts;

    public ConsensusStatistics()
    {
        mDualStrandMismatchReadGroupCount = 0;
        mDualStrandMismatchReadCount = 0;
        mOutcomeCounts = new int[ConsensusOutcome.values().length];
    }

    public void registerOutcome(final ConsensusOutcome outcome)
    {
        ++mOutcomeCounts[outcome.ordinal()];
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

        for(int i = 0; i < mOutcomeCounts.length; ++i)
        {
            mOutcomeCounts[i] += other.mOutcomeCounts[i];
        }
    }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(", ");
        for(ConsensusOutcome outcome : ConsensusOutcome.values())
        {
            if(mOutcomeCounts[outcome.ordinal()] > 0)
                sj.add(format("%s=%d", outcome, mOutcomeCounts[outcome.ordinal()]));
        }

        return format("dualStrandMismatches(groups=%d reads=%s) outcomes(%s)",
                mDualStrandMismatchReadGroupCount, mDualStrandMismatchReadCount, sj);
    }
}