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
    private int mDuplicateQCPassed;
    private int mDuplicateQCFailed;
    private int mPrimaryDuplicateQCPassed;
    private int mPrimaryDuplicateQCFailed;
    private int mMappedQCPassed;
    private int mMappedQCFailed;
    private int mPrimaryMappedQCPassed;
    private int mPrimaryMappedQCFailed;
    private int mPairedQCPassed;
    private int mPairedQCFailed;

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
        mDuplicateQCPassed = 0;
        mDuplicateQCFailed = 0;
        mPrimaryDuplicateQCPassed = 0;
        mPrimaryDuplicateQCFailed = 0;
        mMappedQCPassed = 0;
        mMappedQCFailed = 0;
        mPrimaryMappedQCPassed = 0;
        mPrimaryMappedQCFailed = 0;
        mPairedQCPassed = 0;
        mPairedQCFailed = 0;
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
        mDuplicateQCPassed += other.mDuplicateQCPassed;
        mDuplicateQCFailed += other.mDuplicateQCFailed;
        mPrimaryDuplicateQCPassed += other.mPrimaryDuplicateQCPassed;
        mPrimaryDuplicateQCFailed += other.mPrimaryDuplicateQCFailed;
        mMappedQCPassed += other.mMappedQCPassed;
        mMappedQCFailed += other.mMappedQCFailed;
        mPrimaryMappedQCPassed += other.mPrimaryMappedQCPassed;
        mPrimaryMappedQCFailed += other.mPrimaryMappedQCFailed;
        mPairedQCPassed += other.mPairedQCPassed;
        mPairedQCFailed += other.mPairedQCFailed;
    }

    public void processRead(final SAMRecord read)
    {
        final boolean passesQC = !read.getReadFailsVendorQualityCheckFlag();
        final boolean isSecondary = read.isSecondaryAlignment();
        final boolean isSupp = read.getSupplementaryAlignmentFlag();
        final boolean isDuplicate = read.getDuplicateReadFlag();
        final boolean isMapped = !read.getReadUnmappedFlag();
        final boolean isPaired = read.getReadPairedFlag();

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
                if(isDuplicate)
                {
                    mPrimaryDuplicateQCPassed++;
                }
                if(isMapped)
                {
                    mPrimaryMappedQCPassed++;
                }
            }
            else
            {
                mPrimaryQCFailed++;
                if(isDuplicate)
                {
                    mPrimaryDuplicateQCFailed++;
                }
                if(isMapped)
                {
                    mPrimaryMappedQCFailed++;
                }
            }
        }

        if(isDuplicate)
        {
            if(passesQC)
            {
                mDuplicateQCPassed++;
            }
            else
            {
                mDuplicateQCFailed++;
            }
        }

        if(isMapped)
        {
            if(passesQC)
            {
                mMappedQCPassed++;
            }
            else
            {
                mMappedQCFailed++;
            }
        }

        if(isPaired)
        {
            if(passesQC)
            {
                mPairedQCPassed++;
            }
            else
            {
                mPairedQCFailed++;
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

    public int getDuplicateQCPassed()
    {
        return mDuplicateQCPassed;
    }

    public int getDuplicateQCFailed()
    {
        return mDuplicateQCFailed;
    }

    public int getPrimaryDuplicateQCPassed()
    {
        return mPrimaryDuplicateQCPassed;
    }

    public int getPrimaryDuplicateQCFailed()
    {
        return mPrimaryDuplicateQCFailed;
    }

    public int getMappedQCPassed()
    {
        return mMappedQCPassed;
    }

    public int getMappedQCFailed()
    {
        return mMappedQCFailed;
    }

    public int getPrimaryMappedQCPassed()
    {
        return mPrimaryMappedQCPassed;
    }

    public int getPrimaryMappedQCFailed()
    {
        return mPrimaryMappedQCFailed;
    }

    public int getPairedQCPassed()
    {
        return mPairedQCPassed;
    }

    public int getPairedQCFailed()
    {
        return mPairedQCFailed;
    }
}
