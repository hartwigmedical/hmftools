package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.consumesRefOrUnclippedBases;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.S;

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

    private int mCurrentRefPosition;

    public ReadParseState(final SAMRecord read, boolean isForward)
    {
        Read = read;
        mIsForward = isForward;
        mElementCount = read.getCigar().getCigarElements().size();

        reset();
    }

    public void reset()
    {
        mExhausted = false;

        if(mIsForward)
        {
            mReadIndex = 0;
            mCigarIndex = 0;
        }
        else
        {
            mReadIndex = Read.getReadBases().length - 1;
            mCigarIndex = mElementCount - 1;
        }

        mElementLength = Read.getCigar().getCigarElement(mCigarIndex).getLength();
        mElementType = Read.getCigar().getCigarElement(mCigarIndex).getOperator();
        mElementIndex = 0;

        if(mIsForward)
        {
            mCurrentRefPosition = Read.getAlignmentStart();

            if(mElementType == S)
                mCurrentRefPosition -= mElementLength;
        }
        else
        {
            mCurrentRefPosition = Read.getAlignmentEnd();

            if(mElementType == S)
                mCurrentRefPosition += mElementLength;
        }
    }

    public int readIndex() { return mReadIndex; }
    public byte base() { return mExhausted ? NO_BASE : Read.getReadBases()[mReadIndex]; }
    public byte baseQual() { return mExhausted ? NO_BASE : Read.getBaseQualities()[mReadIndex]; }
    public int refPosition() { return mCurrentRefPosition; }

    public boolean beforeUnclippedPosition(int refPosition)
    {
        return mIsForward ? mCurrentRefPosition > refPosition : mCurrentRefPosition < refPosition;
    }

    public CigarOperator elementType() { return mElementType; }
    public int elementLength() { return mElementLength; }
    public int elementIndex() { return mElementIndex; }
    public int cigarIndex() { return mCigarIndex; }

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
        if(!skipsFirstElement)
        {
            if(mElementType.consumesReadBases())
                mReadIndex += mIsForward ? 1 : -1;

            if(consumesRefOrUnclippedBases(mElementType))
                mCurrentRefPosition += mIsForward ? 1 : -1;
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

    public void moveToRefPosition(int targetPosition)
    {
        if(mIsForward && mCurrentRefPosition >= targetPosition)
            return;
        else if(!mIsForward && mCurrentRefPosition <= targetPosition)
            return;

        while(mCurrentRefPosition != targetPosition && !mExhausted)
        {
            moveNextBase();
        }
    }

    public String toString()
    {
        int effectElementIndex = (mIsForward ? mElementIndex : mElementLength - mElementIndex) + 1;
        return format("index(%d) refPos(%d) cigar(%d: %s element=%d/%d) %s",
                mReadIndex, mCurrentRefPosition, mCigarIndex, mElementType, effectElementIndex, mElementLength,
                mExhausted ? "exhausted" : "active");
    }
}
