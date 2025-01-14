package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class PartitionPerfStats implements Comparable<PartitionPerfStats>
{
    private final String mPartition;
    private double mTotalPc;
    private double mPrePc;
    private double mPostPc;
    private double mConsensusPc;

    public PartitionPerfStats(final String partition, final PerformanceCounter totalPc, final PerformanceCounter prePc,
            final PerformanceCounter postPc, final PerformanceCounter consensusPc)
    {
        mPartition = partition;
        mTotalPc = totalPc.getTotalTime();
        mPrePc = prePc.getTotalTime();
        mPostPc = postPc.getTotalTime();
        mConsensusPc = consensusPc.getTotalTime();
    }

    public PartitionPerfStats(final String partition)
    {
        mPartition = partition;
        mTotalPc = 0.0;
        mPrePc = 0.0;
        mPostPc = 0.0;
        mConsensusPc = 0.0;
    }

    public void merge(final PartitionPerfStats o)
    {
        mTotalPc += o.mTotalPc;
        mPrePc += o.mPrePc;
        mPostPc += o.mPostPc;
        mConsensusPc += o.mConsensusPc;
    }

    @Override
    public String toString()
    {
        return format("PartitionPerfStats(%s) total(%.3f) pre(%.2f%%) pos(%.2f%%) consensus(%.2f%%)", mPartition, mTotalPc, 100.0*mPrePc/mTotalPc, 100.0*mPostPc/mTotalPc, 100.0*mConsensusPc/mTotalPc);
    }

    @Override
    public int compareTo(final PartitionPerfStats o)
    {
        if(mTotalPc > o.mTotalPc)
            return -1;

        if(mTotalPc < o.mTotalPc)
            return 1;

        return 0;
    }
}
