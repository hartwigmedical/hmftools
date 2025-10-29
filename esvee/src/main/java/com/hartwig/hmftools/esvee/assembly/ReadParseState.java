package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.NO_BASE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.DELETE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.INSERT;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.MATCH;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
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

    private double mMismatchPenality;
    public int mHighQualMatches;

    private final List<SequenceDiffInfo> mMismatches;

    public ReadParseState(final boolean moveForward, final Read read, final int startIndex)
    {
        mRead = read;
        mBaseLength = mRead.basesLength();
        mMoveForward = moveForward;
        mStartIndex = startIndex;

        mExhausted = false;
        mReadIndex = 0;
        mIsValid = true;

        // move to the required index
        resetIndex();

        mMismatchPenality = 0;
        mHighQualMatches = 0; // was initialised to 1 for first ref base for ref base building
        mMismatches = Lists.newArrayList();
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
    public List<SequenceDiffInfo> mismatchInfos() { return mMismatches; }
    public void addMismatchInfo(final SequenceDiffInfo seqDiffInfo) { mMismatches.add(seqDiffInfo); }

    public double mismatchPenality() { return mMismatchPenality; }
    public void addMismatch() { ++mMismatchPenality; }

    public void resetMatches()
    {
        mMismatchPenality = 0;
        mHighQualMatches = 0;
        mMismatches.clear();
    }

    public boolean exceedsMaxMismatches(double maxMismatchPenalty) { return mMismatchPenality > maxMismatchPenalty; }

    public void moveOnMatchType(final SequenceDiffInfo seqDiffInfo)
    {
        if(seqDiffInfo.RepeatCount != 0)
            return;

        if(seqDiffInfo.Type == MATCH || seqDiffInfo.Type == BASE)
        {
            moveNext();
        }
        else if(seqDiffInfo.Type == DELETE)
        {
            // skip over base
        }
        else if(seqDiffInfo.Type == INSERT)
        {
            // skip this and the consensus base
            moveNextBases(2);
        }
    }

    public void moveNext() { moveNextBases(1); }

    public void moveNextBases(int baseCount)
    {
        if(mExhausted)
            return;

        for(int i = 0; i < baseCount; ++i)
        {
            if(mMoveForward)
                ++mReadIndex;
            else
                --mReadIndex;

            if(mReadIndex < 0 || mReadIndex >= mBaseLength)
            {
                mExhausted = true;
                break;
            }
        }
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
        return format("%s: index(%d/%d) %s hqMatch(%d) mismatch(%.1f)",
                mRead.id(), mReadIndex, mBaseLength - 1, mExhausted ? "exhausted" : "active", mHighQualMatches, mMismatchPenality);
    }

    @VisibleForTesting
    public void resetAll()
    {
        resetIndex();
        resetMatches();
    }
}
