package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.getUnclippedPosition;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.UNSET;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.bam.SamRecordUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ConsensusState
{
    public final boolean IsForward;
    public final String Chromosome;
    public final Map<String, Object> Attributes;
    public byte[] Bases;
    public byte[] BaseQualities;
    public List<CigarElement> CigarElements;

    public int UnclippedPosStart;
    public int UnclippedPosEnd;
    public int AlignmentStart;
    public int AlignmentEnd;
    public int MapQuality;
    public int NumMutations;

    private ConsensusOutcome mOutcome;
    private final RefGenome mRefGenome;

    private int mCurrentCigarElementLength;
    private CigarOperator mCurrentCigarElementOperator;

    public ConsensusState(final boolean isForward, final String chromosome, final RefGenome refGenome)
    {
        IsForward = isForward;
        Chromosome = chromosome;
        Attributes = Maps.newHashMap();
        mRefGenome = refGenome;
        Bases = null;
        BaseQualities = null;
        CigarElements = Lists.newArrayList();

        UnclippedPosStart = 0;
        UnclippedPosEnd = 0;
        AlignmentStart = 0;
        AlignmentEnd = 0;
        MapQuality = 0;
        NumMutations = 0;

        mCurrentCigarElementLength = 0;
        mCurrentCigarElementOperator = null;

        mOutcome = UNSET;
    }

    public ConsensusOutcome outcome() { return mOutcome; }
    public void setOutcome(final ConsensusOutcome outcome) { mOutcome = outcome; }

    public void setBaseLength(int baseLength)
    {
        Bases = new byte[baseLength];
        BaseQualities = new byte[baseLength];
    }

    public void setBoundaries(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();
        int readEnd = read.getAlignmentEnd();

        int unclippedStart = readStart - CigarUtils.leftSoftClipLength(read);
        int unclippedEnd = readEnd + CigarUtils.rightSoftClipLength(read);

        if(UnclippedPosStart == 0)
        {
            UnclippedPosStart = unclippedStart;
            UnclippedPosEnd = unclippedEnd;
            AlignmentStart = readStart;
            AlignmentEnd = readEnd;
        }
        else
        {
            UnclippedPosStart = min(unclippedStart, UnclippedPosStart);
            UnclippedPosEnd = max(unclippedEnd, UnclippedPosEnd);
            AlignmentStart = min(readStart, AlignmentStart);
            AlignmentEnd = max(readEnd, AlignmentEnd);
        }
    }

    public void setBoundaries(int unclippedStart, int unclippedEnd, int readStart, int readEnd)
    {
        UnclippedPosStart = unclippedStart;
        UnclippedPosEnd = unclippedEnd;
        AlignmentStart = readStart;
        AlignmentEnd = readEnd;
    }

    public static int[] setReadPositionStartOffsets(final List<SAMRecord> reads, final int consensusUnclippedPosition, boolean isStart)
    {
        // convention is to return read's position - consensus position
        int[] positionOffsets = new int[reads.size()];

        for(int i = 0; i < reads.size(); ++i)
        {
            SAMRecord read = reads.get(i);

            int readUnclippedPosition = getUnclippedPosition(read, isStart);
            positionOffsets[i] = readUnclippedPosition - consensusUnclippedPosition;
        }

        return positionOffsets;
    }

    public void addCigarElement(int length, final CigarOperator operator)
    {
        // combine with existing if a match on type
        if(mCurrentCigarElementLength > 0 && mCurrentCigarElementOperator != operator)
            addCurrentCigarElement();

        mCurrentCigarElementLength += length;
        mCurrentCigarElementOperator = operator;
    }

    private void addCurrentCigarElement()
    {
        if(mCurrentCigarElementLength == 0)
            return;

        CigarElement element = new CigarElement(mCurrentCigarElementLength, mCurrentCigarElementOperator);

        if(IsForward)
            CigarElements.add(element);
        else
            CigarElements.add(0, element);

        mCurrentCigarElementLength = 0;
        mCurrentCigarElementOperator = null;
    }

    public void finaliseCigar() { addCurrentCigarElement(); }

    public void setNumMutations()
    {
        NumMutations = 0;
        byte[] refBases = mRefGenome.getRefBases(Chromosome, AlignmentStart, AlignmentEnd);

        if(refBases == null) // abort any attempt to set this property
            return;

        int baseIndex = 0;
        int refBaseIndex = 0;
        for(CigarElement cigarElement : CigarElements)
        {
            CigarOperator cigarOp = cigarElement.getOperator();
            int elemLength = cigarElement.getLength();

            if(isIndelOrMismatch(cigarOp))
            {
                NumMutations += elemLength;
            }
            else if(cigarOp == CigarOperator.M)
            {
                for(int i = 0; i < elemLength; i++)
                {
                    if(refBases[refBaseIndex + i] != Bases[baseIndex + i])
                        ++NumMutations;
                }
            }

            if(cigarOp.consumesReadBases())
                baseIndex += elemLength;

            if(cigarOp.consumesReferenceBases())
                refBaseIndex += elemLength;
        }
    }

    protected static boolean deleteOrSplit(final CigarOperator operator)
    {
        return operator == D || operator == N;
    }

    protected static boolean alignedOrClipped(final CigarOperator operator)
    {
        return operator == M || operator.isClipping();
    }

    protected static boolean consumesRefOrUnclippedBases(final CigarOperator operator)
    {
        return operator.consumesReferenceBases() || operator == S;
    }

    protected static boolean isIndelOrMismatch(CigarOperator cigarOp)
    {
        if(cigarOp == CigarOperator.I)
            return true;

        if(cigarOp == CigarOperator.D)
            return true;

        if(cigarOp == CigarOperator.X)
            return true;

        return false;
    }
}
