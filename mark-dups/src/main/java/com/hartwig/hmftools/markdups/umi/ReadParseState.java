package com.hartwig.hmftools.markdups.umi;

import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.umi.BaseBuilder.NO_BASE;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadParseState
{
    public final SAMRecord Read;

    private boolean mIsForward;
    private int mReadIndex;
    private int mCigarIndex;
    private int mCigarElementIndex;
    private int mElementLength;
    private CigarOperator mElementType;

    private boolean mExhausted;
    private boolean mUseCurrentBase;
    private int mSkipCount;

    public ReadParseState(final SAMRecord read, boolean isForward)
    {
        Read = read;
        mIsForward = isForward;

        if(!read.getReadNegativeStrandFlag())
        {
            mReadIndex = 0;
            mCigarIndex = 0;
        }
        else
        {
            mReadIndex = read.getReadBases().length - 1;
            mCigarIndex = read.getCigar().getCigarElements().size() - 1;
        }

        mElementLength = read.getCigar().getCigarElement(mCigarIndex).getLength();
        mElementType = read.getCigar().getCigarElement(mCigarIndex).getOperator();
        mCigarElementIndex = 0;
        mExhausted = false;
        mSkipCount = 0;
        mUseCurrentBase = true;
    }

    public byte currentBase()
    {
        return mExhausted ? NO_BASE : Read.getReadBases()[mReadIndex];
    }

    public byte currentBaseQual()
    {
        return mExhausted ? NO_BASE : Read.getBaseQualities()[mReadIndex];
    }

    public CigarOperator currentElementType()
    {
        return mElementType;
    }

    public int currentElementLength()
    {
        return mElementLength;
    }

    public boolean useCurrentBase()
    {
        return mUseCurrentBase;
    }

    public boolean exhausted()
    {
        return mExhausted;
    }

    public void moveNext()
    {
        if(mExhausted)
        {
            return;
        }

        if(mIsForward)
        {
            ++mReadIndex;

            if(mReadIndex >= Read.getReadBases().length)
            {
                mExhausted = true;
                return;
            }
        }
        else
        {
            --mReadIndex;

            if(mReadIndex < 0)
            {
                mExhausted = true;
                return;
            }
        }

        ++mCigarElementIndex;

        if(mCigarElementIndex >= mElementLength)
        {
            if(mIsForward)
            {
                ++mCigarIndex;
            }
            else
            {
                --mCigarIndex;
            }

            mElementLength = Read.getCigar().getCigarElement(mCigarIndex).getLength();
            mCigarElementIndex = 0;
        }
    }

    public void skipInsert()
    {
        while(mElementType == I && !exhausted())
        {
            moveNext();
        }
    }

    public void skipPosition()
    {
        if(mElementType == M)
        {

        }

    }

    public void pause(int skipCount)
    {
        mSkipCount = skipCount;
    }

    public String toString()
    {
        return format("index(%d) cigar(%d: %s element=%d/%d) %s",
                mReadIndex, mCigarIndex, mElementType, mCigarElementIndex, mElementLength, mExhausted ? "exhausted" : "active");
    }
}
