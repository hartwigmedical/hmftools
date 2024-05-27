package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.consensus.BaseBuilder.NO_BASE;

import static htsjdk.samtools.CigarOperator.I;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadParseState
{
    public final SAMRecord Read;
    private final boolean mIsForward;
    private final int mElementCount;

    private int mReadIndex;
    private int mCigarIndex; // index of the current element amongst the CIGAR
    private int mElementIndex; // index within the current CIGAR element
    private int mElementLength;
    private CigarOperator mElementType;

    private boolean mExhausted;

    public ReadParseState(final SAMRecord read, boolean isForward)
    {
        Read = read;
        mIsForward = isForward;
        mElementCount = read.getCigar().getCigarElements().size();

        if(mIsForward)
        {
            mReadIndex = 0;
            mCigarIndex = 0;
        }
        else
        {
            mReadIndex = read.getReadBases().length - 1;
            mCigarIndex = mElementCount - 1;
        }

        mElementLength = read.getCigar().getCigarElement(mCigarIndex).getLength();
        mElementType = read.getCigar().getCigarElement(mCigarIndex).getOperator();
        mElementIndex = 0;
        mExhausted = false;
    }

    public byte currentBase() { return mExhausted ? NO_BASE : Read.getReadBases()[mReadIndex]; }
    public byte currentBaseQual() { return mExhausted ? NO_BASE : Read.getBaseQualities()[mReadIndex]; }

    public CigarOperator elementType() { return mElementType; }
    public int elementLength() { return mElementLength; }
    public int elementIndex() { return mElementIndex; }

    public boolean exhausted() { return mExhausted; }

    public void moveNextBase()
    {
        if(mExhausted)
            return;

        ++mElementIndex;

        boolean skipsFirstElement = false;

        if(mElementIndex >= mElementLength)
        {
            skipsFirstElement = !mElementType.consumesReadBases() && (mCigarIndex == 0 || mCigarIndex == mElementCount - 1);

            if(mIsForward)
                ++mCigarIndex;
            else
                --mCigarIndex;

            if(mCigarIndex < 0 || mCigarIndex >= Read.getCigar().getCigarElements().size())
            {
                mExhausted = true;
                return;
            }

            mElementLength = Read.getCigar().getCigarElement(mCigarIndex).getLength();
            mElementType = Read.getCigar().getCigarElement(mCigarIndex).getOperator();
            mElementIndex = 0;
        }

        // move the read index to the start of the new element
        if(!skipsFirstElement && mElementType.consumesReadBases())
        {
            if(mIsForward)
                ++mReadIndex;
            else
                --mReadIndex;
        }
    }

    public void skipInsertElement()
    {
        while(mElementType == I && !exhausted())
        {
            moveNextBase();
        }
    }

    public void skipNonReadBaseElement()
    {
        if(mExhausted)
            return;

        boolean skipsFirstElement = !mElementType.consumesReadBases() && (mCigarIndex == 0 || mCigarIndex == mElementCount - 1);

        if(mIsForward)
            ++mCigarIndex;
        else
            --mCigarIndex;

        if(mCigarIndex < 0 || mCigarIndex >= Read.getCigar().getCigarElements().size())
        {
            mExhausted = true;
            return;
        }

        mElementLength = Read.getCigar().getCigarElement(mCigarIndex).getLength();
        mElementType = Read.getCigar().getCigarElement(mCigarIndex).getOperator();
        mElementIndex = 0;

        // if cigar starts with a H, no need to move to the next base
        if(!skipsFirstElement && mElementType.consumesReadBases())
        {
            if(mIsForward)
                ++mReadIndex;
            else
                --mReadIndex;
        }
    }

    public String toString()
    {
        return format("index(%d) cigar(%d: %s element=%d/%d) %s",
                mReadIndex, mCigarIndex, mElementType, mElementIndex, mElementLength, mExhausted ? "exhausted" : "active");
    }
}
