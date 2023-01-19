package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.UNSET;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class IndelReadUtils
{
    public ConsensusOutcome buildIndelComponents(final List<SAMRecord> reads, boolean isForward, final ConsensusState consensusState)
    {
        List<ReadState> readStates = Lists.newArrayListWithExpectedSize(reads.size());

        int maxIndelCount = 0;

        for(SAMRecord read : reads)
        {
            readStates.add(new ReadState(read));

            int indelCount =
                    read.getCigar().getCigarElements().stream().mapToInt(x -> x.getOperator() == I || x.getOperator() == D ? 1 : 0).sum();

            maxIndelCount = max(maxIndelCount, indelCount);
        }


        int baseLength = consensusState.Bases.length;

        int baseIndex = isForward ? 0 : baseLength - 1;

        ConsensusOutcome outcome = UNSET;

        while(baseIndex >= 0 && baseIndex < baseLength)
        {
            // update the state for each read
            for(ReadState readState : readStates)
            {
                readState.moveNext(isForward);
            }

            // check bases at this index


            if(isForward)
                ++baseIndex;
            else
                --baseIndex;
        }

        return outcome;
    }

    private static final byte INVALID_BASE = -1;

    private class ReadState
    {
        public final SAMRecord Read;

        private int mReadIndex;
        private int mCigarIndex;
        private int mCigarElementIndex;
        private int mCigarElementLength;
        private boolean mExhausted;

        public ReadState(final SAMRecord read)
        {
            Read = read;

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

            mCigarElementLength = read.getCigar().getCigarElement(mCigarIndex).getLength();
            mCigarElementIndex = 0;
            mExhausted = false;
        }

        public boolean exhausted() { return mExhausted; }

        public void moveNext(boolean isForward)
        {
            if(mExhausted)
                return;

            if(isForward)
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

            if(mCigarElementIndex >= mCigarElementLength)
            {
                if(isForward)
                    ++mCigarIndex;
                else
                    --mCigarIndex;

                mCigarElementLength = Read.getCigar().getCigarElement(mCigarIndex).getLength();
                mCigarElementIndex = 0;
            }
        }

        public byte currentBase() { return mExhausted ? INVALID_BASE : Read.getReadBases()[mReadIndex]; }
        public byte currentBaseQual() { return mExhausted ? INVALID_BASE : Read.getBaseQualities()[mReadIndex]; }
        public CigarOperator currentCigarType() { return Read.getCigar().getCigarElement(mCigarIndex).getOperator(); }

        public String toString()
        {
            return format("index(%d) cigar(%d element=%d/%d) %s",
                    mReadIndex, mCigarIndex, mCigarElementIndex, mCigarElementLength, mExhausted ? "exhausted" : "active");
        }
    }

}
