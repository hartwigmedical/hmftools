package com.hartwig.hmftools.bamtools.metrics;

import htsjdk.samtools.SAMRecord;

public class FlagStats
{
    private int mTotalQCPassed;
    private int mTotalQCFailed;
    private int mPrimaryQCPassed;
    private int mPrimaryQCFailed;
    private int mSecondaryQCPassed;
    private int mSecondaryQCFailed;
    private int mSuppQCPassed;
    private int mSuppQCFailed;

    public FlagStats()
    {
        mTotalQCPassed = 0;
        mTotalQCFailed = 0;
        mPrimaryQCPassed = 0;
        mPrimaryQCFailed = 0;
        mSecondaryQCPassed = 0;
        mSecondaryQCFailed = 0;
        mSuppQCPassed = 0;
        mSuppQCFailed = 0;
    }

    public void merge(final FlagStats other)
    {
        mTotalQCPassed += other.mTotalQCPassed;
        mTotalQCFailed += other.mTotalQCFailed;
        mPrimaryQCPassed += other.mPrimaryQCPassed;
        mPrimaryQCFailed += other.mPrimaryQCFailed;
        mSecondaryQCPassed += other.mSecondaryQCPassed;
        mSecondaryQCFailed += other.mSecondaryQCFailed;
        mSuppQCPassed += other.mSuppQCPassed;
        mSuppQCFailed += other.mSuppQCFailed;
    }

    public void processRead(final SAMRecord read)
    {
        boolean passesQC = !read.getReadFailsVendorQualityCheckFlag();
        boolean isSecondary = read.isSecondaryAlignment();
        boolean isSupp = read.getSupplementaryAlignmentFlag();

        if(passesQC)
        {
            mTotalQCPassed++;
        }
        else
        {
            mTotalQCFailed++;
        }

        if(isSecondary)
        {
            if(passesQC)
            {
                mSecondaryQCPassed++;
            }
            else
            {
                mSecondaryQCFailed++;
            }
        }
        else if(isSupp)
        {
            if(passesQC)
            {
                mSuppQCPassed++;
            }
            else
            {
                mSuppQCFailed++;
            }
        }
        else
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

    public int getSecondaryQCPassed()
    {
        return mSecondaryQCPassed;
    }

    public int getSecondaryQCFailed()
    {
        return mSecondaryQCFailed;
    }

    public int getSuppQCPassed()
    {
        return mSuppQCPassed;
    }

    public int getSuppQCFailed()
    {
        return mSuppQCFailed;
    }
}
