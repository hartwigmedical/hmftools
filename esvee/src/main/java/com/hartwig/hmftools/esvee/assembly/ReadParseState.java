package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.NO_BASE;

import com.hartwig.hmftools.esvee.assembly.read.Read;

public class ReadParseState
{
    private final Read mRead;
    private final boolean mMoveForward;
    private final int mBaseLength;
    private final int mStartIndex;
    private int mReadIndex;

    private boolean mExhausted;
    private boolean mIsValid;

    private int mMismatches;
    private int mIndelMismatches;
    public int mHighQualMatches;

    public ReadParseState(final boolean moveForward, final Read read, final int startIndex)
    {
        mRead = read;
        mBaseLength = mRead.basesLength();
        mMoveForward = moveForward;
        mStartIndex = startIndex;

        mExhausted = false;
        mIsValid = true;
        mReadIndex = 0;

        // move to the required index
        resetIndex();

        mMismatches = 0;
        mIndelMismatches = 0;
        mHighQualMatches = 0; // was initialised to 1 for first ref base for ref base building
    }

    public Read read()
    {
        return mRead;
    }

    public boolean isValid() { return mIsValid; }
    public void markInvalid() { mIsValid = false; }

    public boolean exhausted() { return mExhausted; }

    public int readIndex() { return mReadIndex; }

    public byte currentBase() { return mRead.getBases()[mReadIndex]; }
    public byte currentQual() { return mRead.getBaseQuality()[mReadIndex]; }

    public int highQualMatches() { return mHighQualMatches; }
    public void addHighQualMatch() { ++mHighQualMatches; }

    public int mismatches() { return mMismatches; }
    public void addMismatch() { ++mMismatches; }

    public void resetMatches()
    {
        mMismatches = 0;
        mHighQualMatches = 0;
    }

    public int indelMismatches() { return mIndelMismatches; }
    public void addIndelMismatch() { ++mIndelMismatches; }

    public boolean exceedsMaxMismatches(int permittedMismatches) { return mMismatches > permittedMismatches; }

    public void moveNext()
    {
        if(mExhausted)
            return;

        if(mMoveForward)
        {
            ++mReadIndex;
        }
        else
        {
            --mReadIndex;
        }

        mExhausted = mReadIndex < 0 || mReadIndex >= mBaseLength;

        if(mReadIndex < 0 || mReadIndex >= mBaseLength)
            mIsValid = false;
    }

    public void resetIndex()
    {
        mExhausted = false;

        if(mMoveForward)
        {
            mReadIndex = 0;
        }
        else
        {
            mReadIndex = mBaseLength - 1;
        }

        if(mStartIndex >= 0)
        {
            while(mReadIndex != mStartIndex && !mExhausted)
            {
                moveNext();
            }
        }
    }

    public byte getPreviousBase()  { return getAdjacentBase(false); }
    public byte getNextBase() { return getAdjacentBase(true); }

    private byte getAdjacentBase(boolean getUpper)
    {
        if(mMoveForward == getUpper)
            return mReadIndex >= mBaseLength - 1 ? NO_BASE : mRead.getBases()[mReadIndex + 1];
        else
            return mReadIndex < 1 ? NO_BASE : mRead.getBases()[mReadIndex - 1];
    }

    public byte[] getPreviousBases(int baseCount) { return getAdjacentBases(baseCount, !mMoveForward); }
    public byte[] getNextBases(int baseCount) { return getAdjacentBases(baseCount, mMoveForward); }

    private byte[] getAdjacentBases(int baseCount, boolean getUpper)
    {
        byte[] bases = new byte[baseCount];
        int readIndex = mReadIndex + (getUpper ? 1 : -1);

        if(getUpper)
        {
            for(int i = 0; i < baseCount; ++i)
            {
                if(readIndex >= mBaseLength)
                    return null;

                bases[i] = mRead.getBases()[readIndex];
                ++readIndex;
            }
        }
        else
        {
            for(int i = baseCount - 1; i >= 0; --i)
            {
                if(readIndex < 0)
                    return null;

                bases[i] = mRead.getBases()[readIndex];
                --readIndex;
            }
        }

        return bases;
    }

    public String toString()
    {
        return format("%s: index(%d/%d) %s match(%d) mismatch(%d indel=%d)",
                mRead.id(), mReadIndex, mBaseLength, mIsValid ? (mExhausted ? "exhausted" : "active") : "invalid",
                mHighQualMatches, mMismatches, mIndelMismatches);
    }
}
