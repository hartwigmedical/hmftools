package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;

import com.hartwig.hmftools.esvee.assembly.read.Read;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class ReadParseState
{
    private final Read mRead;
    private final int mJunctionIndex;
    private final boolean mMoveForward;

    private int mReadIndex;
    private int mRefPosition;

    private final int mElementCount;
    private final List<CigarElement> mElements;
    private int mCigarIndex; // index of the current element amongst the CIGAR
    private CigarElement mElement;
    private int mElementIndex; // index within the current CIGAR element

    private boolean mExhausted;
    private boolean mIsValid;
    private int mRefBaseCount;

    private final int mPermittedMismatches;
    private int mMismatches;
    private int mIndelMismatches;
    public int mHighQualMatches;

    public ReadParseState(final boolean isForwardJunction, final Read read, final int junctionIndex, final int refBaseCount)
    {
        mRead = read;
        mJunctionIndex = junctionIndex;
        mMoveForward = !isForwardJunction;
        mRefBaseCount = refBaseCount;

        mElements = mRead.cigarElements();
        mElementCount = mElements.size();
        mExhausted = false;
        mIsValid = true;
        mReadIndex = 0;
        mRefPosition = 0;
        mCigarIndex = 0;
        mElement = null;
        mElementIndex = 0;

        // move to the junction index and corresponding cigar element
        resetIndex();

        mPermittedMismatches = mismatchesPerComparisonLength(mRefBaseCount);
        mMismatches = 0;
        mIndelMismatches = 0;
        mHighQualMatches = 1; // assumes junction ref base matches
    }

    public Read read()
    {
        return mRead;
    }

    public boolean isValid() { return mIsValid; }
    public void markInvalid() { mIsValid = false; }

    public boolean exhausted() { return mExhausted; }

    public int refPosition() { return mRefPosition; }
    public int readIndex() { return mReadIndex; }

    public CigarOperator operator() { return mElement.getOperator(); }
    public int elementLength() { return mElement.getLength(); }
    public byte currentBase() { return mRead.getBases()[mReadIndex]; }
    public byte currentQual() { return mRead.getBaseQuality()[mReadIndex]; }

    public int highQualMatches() { return mHighQualMatches; }
    public void addHighQualMatch() { ++mHighQualMatches; }

    public int mismatches() { return mMismatches; }
    public void addMismatch() { ++mMismatches; }

    public int indelMismatches() { return mIndelMismatches; }
    public void addIndelMismatch() { ++mIndelMismatches; }

    public boolean exceedsMaxMismatches() { return mMismatches > mPermittedMismatches; }

    public void moveStart()
    {
        // skips hard-clips if present
        if(mElement.getOperator() != H)
            return;

        if(mMoveForward)
        {
            ++mCigarIndex;

            if(mCigarIndex >= mElementCount)
            {
                mExhausted = true;
                return;
            }

            mElementIndex = 0;
            mElement = mElements.get(mCigarIndex);
        }
        else
        {
            --mCigarIndex;

            if(mCigarIndex < 0)
            {
                mExhausted = true;
                return;
            }

            mElement = mElements.get(mCigarIndex);
            mElementIndex = mElement.getLength() - 1;
        }
    }

    public void moveNext()
    {
        if(mExhausted)
            return;

        if(mMoveForward)
        {
            ++mElementIndex;

            if(mElementIndex >= mElement.getLength())
            {
                // determine whether to move index and position base on previous element
                boolean wasClipping = mElement.getOperator().isClipping();

                ++mCigarIndex;

                if(mCigarIndex >= mElementCount)
                {
                    mExhausted = true;
                    return;
                }

                mElementIndex = 0;
                mElement = mElements.get(mCigarIndex);

                if(mCigarIndex > 0 && mElement.getOperator().isClipping()) // stop at a soft-clip
                {
                    mExhausted = true;
                }
                else
                {
                    if(!wasClipping && mElement.getOperator() != I)
                        ++mRefPosition;

                    if(mElement.getOperator().consumesReadBases())
                        ++mReadIndex;
                }
            }
            else
            {
                if(mElement.getOperator().consumesReferenceBases())
                    ++mRefPosition;

                if(mElement.getOperator().consumesReadBases())
                    ++mReadIndex;
            }
        }
        else
        {
            --mElementIndex;

            if(mElementIndex < 0 || mElement.getOperator() == H)
            {
                boolean wasClipping = mElement.getOperator().isClipping();

                --mCigarIndex;

                if(mCigarIndex < 0)
                {
                    mExhausted = true;
                    return;
                }

                mElement = mElements.get(mCigarIndex);
                mElementIndex = mElement.getLength() - 1;

                if(mCigarIndex < mElementCount - 1 && mElement.getOperator().isClipping())
                {
                    mExhausted = true;
                }
                else
                {
                    if(!wasClipping && mElement.getOperator() != I)
                        --mRefPosition;

                    if(mElement.getOperator().consumesReadBases())
                        --mReadIndex;
                }
            }
            else
            {
                if(mElement.getOperator().consumesReferenceBases())
                    --mRefPosition;

                if(mElement.getOperator().consumesReadBases())
                    --mReadIndex;
            }
        }

        if(mReadIndex < 0 || mReadIndex >= mRead.getBases().length)
            mIsValid = false;
    }

    public void skipInsert()
    {
        while(!exhausted() && operator() == I)
        {
            moveNext();
        }
    }

    public void resetIndex()
    {
        mExhausted = false;

        if(mMoveForward)
        {
            mReadIndex = 0;
            mRefPosition = mRead.alignmentStart();
            mCigarIndex = 0;
            mElement = mElements.get(mCigarIndex);
            mElementIndex = 0;
        }
        else
        {
            mReadIndex = mRead.basesLength() - 1;
            mRefPosition = mRead.alignmentEnd();
            mCigarIndex = mElementCount - 1;
            mElement = mElements.get(mCigarIndex);
            mElementIndex = mElement.getLength() - 1;
        }

        moveStart();

        while(mReadIndex != mJunctionIndex && !mExhausted)
        {
            moveNext();
        }
    }

    public String toString()
    {
        return format("%s: cigar(%s) refBases(%d) index(%d junc=%d) refPos(%d) element(%d/%s %d/%d) %s match(%d) mismatch(%d indel=%d)",
                mRead.id(), mRead.cigarString(), mRefBaseCount, mReadIndex, mJunctionIndex, mRefPosition,
                mElementIndex + 1, mElement.toString(), mCigarIndex + 1, mElementCount,
                mIsValid ? (mExhausted ? "exhausted" : "active") : "invalid",
                mHighQualMatches, mMismatches, mIndelMismatches);
    }
}
