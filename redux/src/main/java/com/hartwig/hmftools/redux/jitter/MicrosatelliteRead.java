package com.hartwig.hmftools.redux.jitter;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.ConsensusType.NONE;
import static com.hartwig.hmftools.redux.ReduxConfig.SEQUENCING_TYPE;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.LOW_BASE_QUAL_FLANKING_BASES;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.MIN_FLANKING_BASE_MATCHES;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class MicrosatelliteRead
{
    private static final ThreadLocal<MicrosatelliteRead> THREAD_INSTANCE = new ThreadLocal<>()
    {
        @Override
        protected MicrosatelliteRead initialValue()
        {
            return new MicrosatelliteRead();
        }
    };

    private ConsensusType mConsensusType;

    private boolean mIsValidRead;

    private int mAlignedBases;
    private int mInsertedBases;
    private int mDeletedBases;

    private int mRepeatUnits;
    private int mRepeatLength;
    private int mJitterLength;

    public MicrosatelliteRead()
    {
        clear();
    }

    public static MicrosatelliteRead from(
            final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record, final ConsensusMarker consensusMarker)
    {
        MicrosatelliteRead instance = THREAD_INSTANCE.get();
        instance.analyse(refGenomeMicrosatellite, record, consensusMarker);
        return instance;
    }

    public boolean isValidRead() { return mIsValidRead; }

    public ConsensusType consensusType() { return mConsensusType; }
    public int readRepeatLength() { return mRepeatLength; }
    public int numRepeatUnits()
    {
        return mRepeatUnits;
    }
    public int jitter() { return mJitterLength; }

    @VisibleForTesting
    public void analyse(final RefGenomeMicrosatellite microsatelliteRepeat, final SAMRecord record, final ConsensusMarker consensusMarker)
    {
        clear();

        // this read needs to wholly contain the repeat to be counted
        int msPosStart = microsatelliteRepeat.referenceStart();
        int msPosEnd = microsatelliteRepeat.referenceEnd();

        // most basic condition is that the read has aligned bases around the repeat including the defined flanks
        if(record.getAlignmentStart() > msPosStart - MIN_FLANKING_BASE_MATCHES
        || record.getAlignmentEnd() < msPosEnd + MIN_FLANKING_BASE_MATCHES)
        {
            return;
        }

        int alignedBeforeRepeat = 0;
        int alignedAfterRepeat = 0;
        int msReadIndexStart = -1;

        // gather key info from the read
        int readIndex = 0;
        int refPosition = record.getAlignmentStart();

        for(int i = 0; i < record.getCigar().getCigarElements().size(); ++i)
        {
            CigarElement element = record.getCigar().getCigarElements().get(i);

            int endRefPos = refPosition + element.getLength() - 1;

            if(msReadIndexStart < 0 && refPosition <= msPosStart && endRefPos >= msPosStart)
                msReadIndexStart = readIndex + msPosStart - refPosition;

            // note soft-clip positions are not checked since the aligned positions vs repeat flanks covers this

            if(element.getOperator() == D)
            {
                if(refPosition >= msPosStart && endRefPos <= msPosEnd)
                {
                    mDeletedBases += element.getLength();
                }
                else if((refPosition < msPosStart && endRefPos >= msPosStart - 1) || (refPosition <= msPosEnd + 1 && endRefPos > msPosEnd))
                {
                    // drop the read if the delete covers the start or the end of the repeat or is just before its start
                    return;
                }
            }
            else if(element.getOperator() == I)
            {
                if(refPosition >= msPosStart && refPosition <= msPosEnd + 1)
                {
                    // check whether repeat unit is a multiple of the repeat
                    for(int j = 0; j < element.getLength(); ++j)
                    {
                        if(record.getReadBases()[readIndex + j] != microsatelliteRepeat.Unit[j % microsatelliteRepeat.Unit.length])
                            return;
                    }

                    mInsertedBases += element.getLength();
                }
                else if(refPosition >= msPosStart - MIN_FLANKING_BASE_MATCHES &&
                        refPosition <= msPosEnd + MIN_FLANKING_BASE_MATCHES + 1)
                {
                    // the insert cannot be within the flanks of the repeat
                    return;
                }
            }
            else if(element.getOperator() == M)
            {
                // capture aligned bases before, across and after the repeat
                int numAlignedToMs = min(endRefPos + 1, msPosEnd + 1) - max(refPosition, msPosStart);

                if(numAlignedToMs > 0)
                {
                    mAlignedBases += numAlignedToMs;
                }

                int numAlignedToFlankingStart = min(endRefPos + 1, msPosStart) - max(refPosition, msPosStart - MIN_FLANKING_BASE_MATCHES);

                if(numAlignedToFlankingStart > 0)
                {
                    alignedBeforeRepeat += numAlignedToFlankingStart;
                }

                int numAlignedToFlankingEnd = min(endRefPos + 1, msPosEnd + 1 + MIN_FLANKING_BASE_MATCHES) - max(refPosition, msPosEnd + 1);

                if(numAlignedToFlankingEnd > 0)
                {
                    alignedAfterRepeat += numAlignedToFlankingEnd;
                }
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();
        }

        if(alignedBeforeRepeat < MIN_FLANKING_BASE_MATCHES || alignedAfterRepeat < MIN_FLANKING_BASE_MATCHES)
            return;

        if(msReadIndexStart < 0)
            return;

        int msReadIndexEnd = msReadIndexStart + msPosEnd - msPosStart;

        if(!hasValidBaseQualities(record, msReadIndexStart, msReadIndexEnd))
            return;

        mIsValidRead = true;

        mConsensusType = consensusMarker.consensusType(microsatelliteRepeat, record);

        mRepeatLength = mAlignedBases + mInsertedBases;
        mRepeatUnits = mRepeatLength / microsatelliteRepeat.Unit.length;
        mJitterLength = mRepeatUnits - microsatelliteRepeat.RepeatCount;
    }

    private void clear()
    {
        mConsensusType = NONE;
        mIsValidRead = false;

        mAlignedBases = 0;
        mInsertedBases = 0;
        mDeletedBases = 0;
        mRepeatUnits = 0;
        mRepeatLength = 0;
        mJitterLength = 0;
    }

    private boolean hasValidBaseQualities(final SAMRecord record, int msReadIndexStart, int msReadIndexEnd)
    {
        int msFlankReadIndexStart = msReadIndexStart - LOW_BASE_QUAL_FLANKING_BASES;
        int msFlankReadIndexEnd = msReadIndexEnd + LOW_BASE_QUAL_FLANKING_BASES;

        if(msFlankReadIndexStart < 0 || msFlankReadIndexEnd >= record.getBaseQualities().length)
            return false;

        for(int i = msFlankReadIndexStart; i <= msFlankReadIndexEnd; ++i)
        {
            if(BaseQualAdjustment.isUncertainBaseQual(record.getBaseQualities()[i]))
                return false;

            if(BaseQualAdjustment.isMediumBaseQual(record.getBaseQualities()[i], SEQUENCING_TYPE))
                return false;
        }

        return true;
    }
}
