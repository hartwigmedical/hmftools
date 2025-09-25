package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.UNSET;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ConsensusState
{
    public final boolean IsForward;
    public final String Chromosome;
    public final Map<String, Object> Attributes;
    private final RefGenome mRefGenome;
    public byte[] Bases;
    public byte[] BaseQualities;
    public List<CigarElement> CigarElements;

    private int mCurrentCigarElementLength;
    private CigarOperator mCurrentCigarElementOperator;

    public int MinUnclippedPosStart;
    public int MaxUnclippedPosEnd;
    public int MinAlignedPosStart;
    public int MaxAlignedPosEnd;
    public int MapQuality;
    public int NumMutations;

    private ConsensusOutcome mOutcome;

    public ConsensusState(final boolean isForward, final String chromosome, final RefGenome refGenome)
    {
        IsForward = isForward;
        Chromosome = chromosome;
        Attributes = Maps.newHashMap();
        mRefGenome = refGenome;
        Bases = null;
        BaseQualities = null;
        CigarElements = Lists.newArrayList();

        MinUnclippedPosStart = 0;
        MaxUnclippedPosEnd = 0;
        MinAlignedPosStart = 0;
        MaxAlignedPosEnd = 0;
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
        int unclippedStart = read.getCigar().isLeftClipped() ? readStart - read.getCigar().getFirstCigarElement().getLength() : readStart;
        int unclippedEnd = read.getCigar().isRightClipped() ? readEnd + read.getCigar().getLastCigarElement().getLength() : readEnd;

        if(MinUnclippedPosStart == 0)
        {
            MinUnclippedPosStart = unclippedStart;
            MaxUnclippedPosEnd = unclippedEnd;
            MinAlignedPosStart = readStart;
            MaxAlignedPosEnd = readEnd;
        }
        else
        {
            MinUnclippedPosStart = min(unclippedStart, MinUnclippedPosStart);
            MaxUnclippedPosEnd = max(unclippedEnd, MaxUnclippedPosEnd);
            MinAlignedPosStart = min(readStart, MinAlignedPosStart);
            MaxAlignedPosEnd = max(readEnd, MaxAlignedPosEnd);
        }
    }

    public void setBoundaries(int unclippedStart, int unclippedEnd, int readStart, int readEnd)
    {
        MinUnclippedPosStart = unclippedStart;
        MaxUnclippedPosEnd = unclippedEnd;
        MinAlignedPosStart = readStart;
        MaxAlignedPosEnd = readEnd;
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
        byte[] refBases = mRefGenome.getRefBases(Chromosome, MinAlignedPosStart, MaxAlignedPosEnd);

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

    private static boolean isIndelOrMismatch(CigarOperator cigarOp)
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
