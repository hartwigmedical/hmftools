package com.hartwig.hmftools.bamtools.metrics;

import htsjdk.samtools.SAMRecord;

public class FlagStats
{
    private int mTotalQCPassed;
    private int mTotalQCFailed;
    private int mPrimaryQCPassed;
    private int mPrimaryQCFailed;

    public FlagStats()
    {
        mTotalQCPassed = 0;
        mTotalQCFailed = 0;
        mPrimaryQCPassed = 0;
        mPrimaryQCFailed = 0;
    }

    public void merge(final FlagStats other)
    {
        mTotalQCPassed += other.mTotalQCPassed;
        mTotalQCFailed += other.mTotalQCFailed;
        mPrimaryQCPassed += other.mPrimaryQCPassed;
        mPrimaryQCFailed += other.mPrimaryQCFailed;
    }

    public void processRead(final SAMRecord read)
    {
        boolean passesQC = !read.getReadFailsVendorQualityCheckFlag();

        if(passesQC)
        {
            mTotalQCPassed++;
        }
        else
        {
            mTotalQCFailed++;
        }

        if (!read.isSecondaryOrSupplementary())
        {
            if(passesQC)
            {
                mPrimaryQCPassed++;
            }
            else
            {
                mPrimaryQCFailed++;
            }
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

    public int getPrimaryQCPassed()
    {
        return mPrimaryQCPassed;
    }

    public int getPrimaryQCFailed()
    {
        return mPrimaryQCFailed;
    }
}
