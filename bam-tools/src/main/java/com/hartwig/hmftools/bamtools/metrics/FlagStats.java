package com.hartwig.hmftools.bamtools.metrics;

import htsjdk.samtools.SAMRecord;

public class FlagStats
{
    private int mTotalQCPassed;
    private int mTotalQCFailed;

    public FlagStats()
    {
        mTotalQCPassed = 0;
        mTotalQCFailed = 0;
    }

    public void merge(final FlagStats other)
    {
        mTotalQCPassed += other.mTotalQCPassed;
        mTotalQCFailed += other.mTotalQCFailed;
    }

    public void processRead(final SAMRecord read)
    {
        if(read.getReadFailsVendorQualityCheckFlag())
        {
            mTotalQCFailed++;
        }
        else
        {
            mTotalQCPassed++;
        }
    }

    public int getTotalQCPassed()
    {
        return mTotalQCPassed;
    }

    public int getTotalQCFailed()
    {
        return mTotalQCFailed;
    }
}
