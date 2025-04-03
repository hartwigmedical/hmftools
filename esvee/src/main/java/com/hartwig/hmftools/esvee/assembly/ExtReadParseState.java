package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findLineExtensionEndIndex;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;

import com.hartwig.hmftools.esvee.assembly.read.Read;

public class ExtReadParseState
{
    private final Read mRead;
    private final boolean mIsForwardJunction;
    private final int mJunctionIndex;
    private final int mExtensionLength;
    private int mCurrentIndex;
    private boolean mExhausted;
    private final int mPermittedMismatches;

    private int mMismatches;
    private int mHighQualMatches;
    private final int mLineExtensionIndex;

    public ExtReadParseState(final Read read, final int junctionIndex, final int extensionLength, final boolean isForwardJunction)
    {
        mRead = read;
        mIsForwardJunction = isForwardJunction;
        mJunctionIndex = junctionIndex;
        mCurrentIndex = junctionIndex;
        mExtensionLength = extensionLength;
        mPermittedMismatches = mismatchesPerComparisonLength(extensionLength);
        mExhausted = false;

        mMismatches = 0;
        mHighQualMatches = 0;

        if(read.hasLineTail())
        {
            byte lineBase = mIsForwardJunction ? LINE_BASE_T : LINE_BASE_A;
            int extensionIndex = mIsForwardJunction ? mJunctionIndex + 1 : mJunctionIndex - 1;
            mLineExtensionIndex = findLineExtensionEndIndex(read, lineBase, extensionIndex, mIsForwardJunction);
        }
        else
        {
            mLineExtensionIndex = -1;
        }
    }

    public Read read()
    {
        return mRead;
    }

    public int junctionIndex()
    {
        return mJunctionIndex;
    }
    public int extensionLength()
    {
        return mExtensionLength;
    }

    public boolean exhausted() { return mExhausted; }

    public void resetMatches()
    {
        mMismatches = 0;
        mHighQualMatches = 0;
    }

    public byte currentBase() { return mRead.getBases()[mCurrentIndex]; }
    public byte currentQual() { return mRead.getBaseQuality()[mCurrentIndex]; }
    public int currentIndex() { return mCurrentIndex; }

    public void moveNext()
    {
        if(mIsForwardJunction)
        {
            ++mCurrentIndex;
        }
        else
        {
            --mCurrentIndex;
        }

        mExhausted = mCurrentIndex < 0 || mCurrentIndex >= mRead.basesLength();
    }

    public void resetIndex()
    {
        mCurrentIndex = mJunctionIndex;
        mExhausted = false;
    }

    public void skipBases(int skipBaseCount)
    {
        // skip ahead of this read's extra repeats
        for(int i = 0; i < skipBaseCount; ++i)
        {
            moveNext();

            if(exhausted())
                break;
        }
    }

    public int matchedBases() { return mExtensionLength - mMismatches; }
    public int mismatches() { return mMismatches; }
    public int highQualMatches() { return mHighQualMatches; }

    public void addMismatch() { ++mMismatches; }

    public void addHighQualMatch()
    {
        if(mCurrentIndex != mJunctionIndex)
        {
            ++mHighQualMatches;
        }
    }

    public int permittedMismatches() { return mPermittedMismatches; }
    public boolean exceedsMaxMismatches() { return mMismatches > mPermittedMismatches; }

    public void movePastLineExtension(byte lineBase, boolean countMismatches)
    {
        if(mLineExtensionIndex < 0)
        {
            return;
        }

        while(!exhausted() && mCurrentIndex != mLineExtensionIndex)
        {
            byte base = currentBase();
            boolean aboveMinQual = aboveMinQual(currentQual());

            if(base != lineBase)
            {
                if(aboveMinQual)
                {
                    if(countMismatches)
                    {
                        ++mMismatches;
                    }
                }
            }
            else
            {
                if(countMismatches && aboveMinQual)
                {
                    ++mHighQualMatches;
                }
            }

            moveNext();
        }

        // move 1 more base
        moveNext();
    }

    public int countRefRepeat(final byte[] repeatBytes)
    {
        int juncIndex = mJunctionIndex;
        int repeatLength = repeatBytes.length;
        int repeatCount = 0;

        if(mIsForwardJunction)
        {
            while(true)
            {
                boolean matched = true;
                int index = juncIndex - repeatLength + 1;
                for(int i = 0; i < repeatLength; ++i)
                {
                    if(index + i < 0)
                    {
                        matched = false;
                        break;
                    }
                    else if(mRead.getBases()[index + i] != repeatBytes[i])
                    {
                        matched = false;
                        break;
                    }
                }

                if(!matched)
                    break;

                ++repeatCount;
                juncIndex -= repeatLength;
            }
        }
        else
        {
            while(true)
            {
                boolean matched = true;
                int index = juncIndex;
                for(int i = 0; i < repeatLength; ++i)
                {
                    if(index + i >= mRead.getBases().length)
                    {
                        matched = false;
                        break;
                    }
                    else if(mRead.getBases()[index + i] != repeatBytes[i])
                    {
                        matched = false;
                        break;
                    }
                }

                if(!matched)
                    break;

                ++repeatCount;
                juncIndex += repeatLength;
            }
        }

        return repeatCount;
    }

    public String toString()
    {
        return format("%s: range(%d - %d) cigar(%s) extLen(%d) index(junc=%d cur=%d) match(hq=%d mismatch %d/%d)%s%s",
                mRead.id(), mRead.unclippedStart(), mRead.unclippedEnd(), mRead.cigarString(),
                mExtensionLength, mJunctionIndex, mCurrentIndex, mHighQualMatches, mMismatches, mPermittedMismatches,
                mLineExtensionIndex >= 0 ? format(" lineIndex(%d)", mLineExtensionIndex) : "", mExhausted ? " exhausted" : "");
    }

    // unused - may use something like this for other seq-tech
    public void movePastLowQualBases()
    {
        moveNext();

        while(!exhausted() && belowMinQual(currentQual()))
        {
            moveNext();
        }
    }

    public byte nextHighQualBase()
    {
        int index = mCurrentIndex;
        while(!exhausted() && belowMinQual(mRead.getBaseQuality()[index]))
        {
            ++index;
        }

        if(exhausted())
            return -1;
        else
            return mRead.getBases()[index];
    }
}
