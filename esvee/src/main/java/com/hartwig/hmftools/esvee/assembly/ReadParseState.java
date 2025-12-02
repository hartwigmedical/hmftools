package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.isLowBaseQual;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.NO_BASE;
import static com.hartwig.hmftools.esvee.assembly.SequenceBuilder.calcMismatchPenalty;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.DELETE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.INSERT;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.MATCH;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.REPEAT;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isHighBaseQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isMediumBaseQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.isUltima;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.I;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class ReadParseState
{
    private final Read mRead;
    private final boolean mMoveForward;
    private final int mBaseLength;
    private final int mStartIndex;
    private int mReadIndex;

    private boolean mExhausted;
    private boolean mValid;

    // ref and cigar position tracking - optional
    private int mRefPosition;
    private final int mElementCount;
    private final List<CigarElement> mElements;
    private int mCigarIndex; // index of the current element amongst the CIGAR
    private CigarElement mElement;
    private int mElementIndex; // index within the current CIGAR element

    // track agreement and mismatches vs consensus sequence
    private boolean mMismatched;
    public int mBaseMatches;
    public int mHighQualMatches; // includes medium qual
    private List<SequenceDiffInfo> mMismatches;

    private List<Integer> mLowQualIndices; // only used for Ultima for now

    public ReadParseState(final boolean moveForward, final Read read, final int startIndex)
    {
        this(moveForward, read, startIndex, false);
    }

    public ReadParseState(final boolean moveForward, final Read read, final int startIndex, boolean trackCigar)
    {
        mRead = read;
        mBaseLength = mRead.basesLength();
        mMoveForward = moveForward;
        mStartIndex = startIndex;

        mExhausted = false;
        mValid = true;
        mReadIndex = 0;

        mRefPosition = 0;
        mCigarIndex = 0;
        mElement = null;
        mElementIndex = 0;

        if(trackCigar)
        {
            mElements = mRead.cigarElements();
            mElementCount = mElements.size();
        }
        else
        {
            mElements = null;
            mElementCount = 0;
        }

        // move to the required index
        resetIndex();

        mMismatched = false;
        mBaseMatches = 0;
        mHighQualMatches = 0; // was initialised to 1 for first ref base for ref base building
        mMismatches = null;
        mLowQualIndices = null;
    }

    public Read read()
    {
        return mRead;
    }

    public boolean mismatched() { return mMismatched; }
    public void markMismatched() { mMismatched = true; }

    public boolean isValid() { return mValid; }
    public void markInvalid() { mValid = true; }

    public boolean exhausted() { return mExhausted; }

    public int startIndex() { return mStartIndex; }
    public int readIndex() { return mReadIndex; }

    public byte currentBase() { return mRead.getBases()[mReadIndex]; }
    public byte currentQual() { return mRead.getBaseQuality()[mReadIndex]; }

    public int refPosition() { return mRefPosition; }
    public CigarOperator operator() { return mElement != null ? mElement.getOperator() : null; }
    public int elementLength() { return mElement != null ? mElement.getLength() : 0; }

    public int overlapBaseCount() { return mMoveForward ? mBaseLength - mStartIndex : mStartIndex + 1; }
    public int evaluatedBaseCount() { return (mMoveForward ? mReadIndex - mStartIndex : mStartIndex - mReadIndex) + 1; }

    public void moveOnMatchType(final SequenceDiffInfo seqDiffInfo)
    {
        // handle differences other than repeat adjustments
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
            moveNextBases(2); // CHECK: would need to know insert (diff) length if novel INDEls < 1 base are supported
        }
    }

    public void moveNextBases(int baseCount)
    {
        if(mExhausted)
            return;

        for(int i = 0; i < baseCount; ++i)
        {
            moveNext();
        }
    }

    public void moveNext()
    {
        if(mExhausted)
            return;

        if(mElements != null)
            moveNextCigarAware();
        else
            moveNextReadIndex();
    }

    private void moveNextReadIndex()
    {
        if(mMoveForward)
            ++mReadIndex;
        else
            --mReadIndex;

        if(mReadIndex < 0 || mReadIndex >= mBaseLength)
        {
            mExhausted = true;
        }
    }

    private void moveNextCigarAware()
    {
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
            mExhausted = true;
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
        }
        else
        {
            mReadIndex = mBaseLength - 1;
        }

        if(mElements != null)
        {
            if(mMoveForward)
            {
                mRefPosition = mRead.alignmentStart();
                mCigarIndex = 0;
                mElement = mElements.get(mCigarIndex);
                mElementIndex = 0;
            }
            else
            {
                mRefPosition = mRead.alignmentEnd();
                mCigarIndex = mElementCount - 1;
                mElement = mElements.get(mCigarIndex);
                mElementIndex = mElement.getLength() - 1;
            }

            checkSkipHardClips();
        }

        if(mStartIndex >= 0)
        {
            while(mReadIndex != mStartIndex && !mExhausted)
            {
                moveNext();
            }
        }
    }

    private byte getAdjacentBase(boolean getUpper)
    {
        if(mMoveForward == getUpper)
            return mReadIndex >= mBaseLength - 1 ? NO_BASE : mRead.getBases()[mReadIndex + 1];
        else
            return mReadIndex < 1 ? NO_BASE : mRead.getBases()[mReadIndex - 1];
    }

    public byte[] getNextBases(int baseCount) { return getAdjacentBases(baseCount, mMoveForward); }
    public byte[] getPreviousBases(int baseCount) { return getAdjacentBases(baseCount, !mMoveForward); }
    public byte getPreviousBase()  { return getAdjacentBase(false); }
    public byte getNextBase() { return getAdjacentBase(true); }

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

    public int nonLowQualBaseCount(int readIndexStart, int readIndexEnd)
    {
        int nonLowQualBases = 0;
        for(int i = readIndexStart; i <= readIndexEnd; ++i)
        {
            if(i >= 0 && i < mRead.getBaseQuality().length)
            {
                if(aboveMinQual(mRead.getBaseQuality()[i]))
                    ++nonLowQualBases;
            }
        }

        return nonLowQualBases;
    }

    public BaseQualType qualType(int index) { return rangeMinQualType(index, index); }

    public BaseQualType straddlingQualType(int readIndexStart, int readIndexEnd)
    {
        BaseQualType qualTypeStart = qualType(readIndexStart);
        BaseQualType qualTypeEnd = qualType(readIndexEnd);

        if(qualTypeStart == qualTypeEnd)
            return qualTypeStart;

        // returns the lower qual-type value
        return qualTypeStart.ordinal() > qualTypeEnd.ordinal() ? qualTypeStart : qualTypeEnd;
    }

    public BaseQualType rangeMinQualType(int readIndexStart, int readIndexEnd)
    {
        if(readIndexStart < 0)
            return BaseQualType.LOW;

        // take the lowest qual across a range
        if(isUltima())
        {
            if(mLowQualIndices == null)
                mLowQualIndices = UltimaBamUtils.extractLowQualIndices(mRead.bamRecord());

            for(int i = readIndexStart; i <= readIndexEnd; ++i)
            {
                if(mLowQualIndices.contains(i))
                    return BaseQualType.LOW;
            }

            return BaseQualType.HIGH;
        }
        else
        {
            final byte[] baseQuals = mRead.getBaseQuality();

            if(readIndexEnd >= baseQuals.length)
                return BaseQualType.LOW;

            boolean hasMedium = false;

            for(int i = readIndexStart; i <= readIndexEnd; ++i)
            {
                if(isLowBaseQual(baseQuals[i]))
                    return BaseQualType.LOW;

                hasMedium |= isMediumBaseQual(baseQuals[i]);
            }

            return hasMedium ? BaseQualType.MEDIUM : BaseQualType.HIGH;
        }
    }

    public int matchedBases() { return mBaseMatches; }

    public void addBaseMatch(boolean isHighQual)
    {
        ++mBaseMatches;

        if(isHighQual)
            ++mHighQualMatches;
    }

    public void addBaseMatches(int count, int highQualCount)
    {
        mBaseMatches += count;
        mHighQualMatches += highQualCount;
    }

    public int highQualMatches() { return mHighQualMatches; }

    public List<SequenceDiffInfo> mismatches() { return mMismatches != null ? mMismatches : Collections.emptyList(); }

    public int mismatchCount(boolean excludeLowQual)
    {
        if(mMismatches == null)
            return 0;

        if(excludeLowQual)
            return (int)mMismatches.stream().filter(x -> x.MismatchPenalty != 0).count();

        return mMismatches.size();
    }

    public void addMismatchInfo(final SequenceDiffInfo seqDiffInfo)
    {
        if(mMismatches == null)
            mMismatches = Lists.newArrayList();

        mMismatches.add(seqDiffInfo);
    }

    public void addMismatch(int consensusIndex, final SequenceDiffType type)
    {
        SequenceDiffInfo seqDiffInfo = type == BASE ?
                SequenceDiffInfo.fromSnv(this, consensusIndex) :
                new SequenceDiffInfo( readIndex(), consensusIndex, "", type, qualType(readIndex()));

        seqDiffInfo.MismatchPenalty = calcMismatchPenalty(seqDiffInfo, null);
        addMismatchInfo(seqDiffInfo);
    }

    public void resetMatches()
    {
        mHighQualMatches = 0;
        mBaseMatches = 0;

        if(mMismatches != null)
            mMismatches.clear();
    }

    public double mismatchPenalty()
    {
        if(mMismatches == null)
            return 0;

        return mMismatches.stream().mapToDouble(x -> x.MismatchPenalty).sum();
    }

    public String toString()
    {
        String cigarInfo = " ";

        if(mElements != null)
        {
            cigarInfo = format(" refPos(%d) element(%d/%s %d/%d) ",
                    mRefPosition, mElementIndex + 1, mElement.toString(), mCigarIndex + 1, mElementCount);
        }

        return format("%s: index(%d/%d)%sstate(%s) match(%d hq=%d pen=%.1f mm=%d)",
                mRead.id(), mReadIndex, mBaseLength - 1, cigarInfo != null ? cigarInfo : " ",
                mMismatched ? "mismatched" : (mExhausted ? "exhausted" : "active"),
                mBaseMatches, mHighQualMatches, mismatchPenalty(), mMismatches != null ? mMismatches.size() : 0);
    }

    public String mismatchInfo()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        if(mMismatches == null || mMismatches.isEmpty())
            return "EXACT";

        if(mMismatched)
            return "MISMATCHED";

        sj.add(format("HQ=%d", mHighQualMatches));

        int lowQualMismatches = 0;
        int snvs = 0;
        int indels = 0;
        int homopolymers = 0;
        int otherRepeats = 0;

        for(SequenceDiffInfo mismatch : mismatches())
        {
            if(mismatch.MismatchPenalty > 0)
            {
                if(mismatch.Type == BASE)
                {
                    ++snvs;
                }
                else if(mismatch.Type == REPEAT)
                {
                    if(mismatch.Bases.length() == 1)
                        ++homopolymers;
                    else
                        ++otherRepeats;
                }
                else if(mismatch.Type == INSERT || (mismatch.Type == DELETE))
                {
                    ++indels;
                }
            }
            else
            {
                ++lowQualMismatches;
            }
        }

        if(lowQualMismatches > 0)
            sj.add(format("LQ=%d", lowQualMismatches));

        if(snvs > 0)
            sj.add(format("SNV=%d", snvs));

        if(homopolymers > 0)
            sj.add(format("HP=%d", homopolymers));

        if(otherRepeats > 0)
            sj.add(format("OR=%d", otherRepeats));

        if(indels > 0)
            sj.add(format("ID=%d", indels));

        return sj.toString();
    }

    @VisibleForTesting
    public void resetAll()
    {
        resetIndex();
        resetMatches();
    }

    public void checkSkipHardClips()
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
}
