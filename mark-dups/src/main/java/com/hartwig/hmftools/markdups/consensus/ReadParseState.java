package com.hartwig.hmftools.markdups.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.consensus.BaseBuilder.NO_BASE;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.N;

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
    private boolean mElementChanged; // on the last move

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
        mElementChanged = false;
    }

    public byte currentBase() { return mExhausted ? NO_BASE : Read.getReadBases()[mReadIndex]; }
    public byte currentBaseQual() { return mExhausted ? NO_BASE : Read.getBaseQualities()[mReadIndex]; }

    public CigarOperator elementType() { return mElementType; }
    public int elementLength() { return mElementLength; }
    public int elementIndex() { return mElementIndex; }

    public boolean exhausted() { return mExhausted; }

    public void moveNext()
    {
        if(mExhausted)
            return;

        ++mElementIndex;

        if(mElementIndex >= mElementLength)
        {
            mElementChanged = true;

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
        else
        {
            mElementChanged = false;
        }

        boolean moveReadIndex = (mElementType != D && mElementType != N);

        if(moveReadIndex)
        {
            if(mIsForward)
                ++mReadIndex;
            else
                --mReadIndex;
        }
    }

    public void skipInsert()
    {
        while(mElementType == I && !exhausted())
        {
            moveNext();
        }
    }

    public String toString()
    {
        return format("index(%d) cigar(%d: %s element=%d/%d) %s",
                mReadIndex, mCigarIndex, mElementType, mElementIndex, mElementLength, mExhausted ? "exhausted" : "active");
    }
}
