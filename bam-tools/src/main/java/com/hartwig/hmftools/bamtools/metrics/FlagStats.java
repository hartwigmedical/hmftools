package com.hartwig.hmftools.bamtools.metrics;

import htsjdk.samtools.SAMRecord;

public class FlagStats
{
    private final FlagQCStats mTotal;
    private final FlagQCStats mPrimary;
    private final FlagQCStats mSecondary;
    private final FlagQCStats mSupp;
    private final FlagQCStats mDuplicate;
    private final FlagQCStats mPrimaryDuplicate;
    private final FlagQCStats mMapped;
    private final FlagQCStats mPrimaryMapped;
    private final FlagQCStats mPaired;
    private final FlagQCStats mRead1;
    private final FlagQCStats mRead2;
    private final FlagQCStats mProperlyPaired;
    private final FlagQCStats mPairMapped;
    private final FlagQCStats mSingleton;
    private final FlagQCStats mInterChrPairMapped;
    private final FlagQCStats mInterChrPairMapQGE5;

    public FlagStats()
    {
        mTotal = new FlagQCStats();
        mPrimary = new FlagQCStats();
        mSecondary = new FlagQCStats();
        mSupp = new FlagQCStats();
        mDuplicate = new FlagQCStats();
        mPrimaryDuplicate = new FlagQCStats();
        mMapped = new FlagQCStats();
        mPrimaryMapped = new FlagQCStats();
        mPaired = new FlagQCStats();
        mRead1 = new FlagQCStats();
        mRead2 = new FlagQCStats();
        mProperlyPaired = new FlagQCStats();
        mPairMapped = new FlagQCStats();
        mSingleton = new FlagQCStats();
        mInterChrPairMapped = new FlagQCStats();
        mInterChrPairMapQGE5 = new FlagQCStats();
    }

    public void merge(final FlagStats other)
    {
        mTotal.merge(other.mTotal);
        mPrimary.merge(other.mPrimary);
        mSecondary.merge(other.mSecondary);
        mSupp.merge(other.mSupp);
        mDuplicate.merge(other.mDuplicate);
        mPrimaryDuplicate.merge(other.mPrimaryDuplicate);
        mMapped.merge(other.mMapped);
        mPrimaryMapped.merge(other.mPrimaryMapped);
        mPaired.merge(other.mPaired);
        mRead1.merge(other.mRead1);
        mRead2.merge(other.mRead2);
        mProperlyPaired.merge(other.mProperlyPaired);
        mPairMapped.merge(other.mPairMapped);
        mSingleton.merge(other.mSingleton);
        mInterChrPairMapped.merge(other.mInterChrPairMapped);
        mInterChrPairMapQGE5.merge(other.mInterChrPairMapQGE5);
    }

    public void processRead(final SAMRecord read)
    {
        final boolean passesQC = !read.getReadFailsVendorQualityCheckFlag();
        final boolean isSecondary = read.isSecondaryAlignment();
        final boolean isSupp = read.getSupplementaryAlignmentFlag();
        final boolean isDuplicate = read.getDuplicateReadFlag();
        final boolean isMapped = !read.getReadUnmappedFlag();
        final boolean isPaired = read.getReadPairedFlag();
        final boolean isProperPair = read.getProperPairFlag();

        mTotal.record(passesQC);

        if(isDuplicate)
        {
            mDuplicate.record(passesQC);
        }

        if(isMapped)
        {
            mMapped.record(passesQC);
        }

        if(isSecondary)
        {
            mSecondary.record(passesQC);
            return;
        }

        if(isSupp)
        {
            mSupp.record(passesQC);
            return;
        }

        // It is a primary read.
        mPrimary.record(passesQC);
        if(isDuplicate)
        {
            mPrimaryDuplicate.record(passesQC);
        }

        if(isMapped)
        {
            mPrimaryMapped.record(passesQC);
        }

        if(!isPaired)
        {
            return;
        }

        // It is a paired primary read.
        mPaired.record(passesQC);

        if(read.getFirstOfPairFlag())
        {
            mRead1.record(passesQC);
        }

        if(read.getSecondOfPairFlag())
        {
            mRead2.record(passesQC);
        }

        if(!isMapped)
        {
            return;
        }

        // It is a mapped paired primary read.
        if(isProperPair)
        {
            mProperlyPaired.record(passesQC);
        }

        if(read.getMateUnmappedFlag())
        {
            mSingleton.record(passesQC);
            return;
        }

        // It is a mapped paired primary read with a mapped mate.
        mPairMapped.record(passesQC);
        if(read.getReferenceName().equals(read.getMateReferenceName()))
        {
            return;
        }

        // It is a mapped paired primary read with a mapped mate toa different chromosome.
        mInterChrPairMapped.record(passesQC);
        if(read.getMappingQuality() >= 5)
        {
            mInterChrPairMapQGE5.record(passesQC);
        }
    }

    public FlagQCStats getTotal()
    {
        return mTotal;
    }

    public FlagQCStats getPrimary()
    {
        return mPrimary;
    }

    public FlagQCStats getSecondary()
    {
        return mSecondary;
    }

    public FlagQCStats getSupp()
    {
        return mSupp;
    }

    public FlagQCStats getDuplicate()
    {
        return mDuplicate;
    }

    public FlagQCStats getPrimaryDuplicate()
    {
        return mPrimaryDuplicate;
    }

    public FlagQCStats getMapped()
    {
        return mMapped;
    }

    public FlagQCStats getPrimaryMapped()
    {
        return mPrimaryMapped;
    }

    public FlagQCStats getPaired()
    {
        return mPaired;
    }

    public FlagQCStats getRead1()
    {
        return mRead1;
    }

    public FlagQCStats getRead2()
    {
        return mRead2;
    }

    public FlagQCStats getProperlyPaired()
    {
        return mProperlyPaired;
    }

    public FlagQCStats getPairMapped()
    {
        return mPairMapped;
    }

    public FlagQCStats getSingleton()
    {
        return mSingleton;
    }

    public FlagQCStats getInterChrPairMapped()
    {
        return mInterChrPairMapped;
    }

    public FlagQCStats getInterChrPairMapQGE5()
    {
        return mInterChrPairMapQGE5;
    }
}
