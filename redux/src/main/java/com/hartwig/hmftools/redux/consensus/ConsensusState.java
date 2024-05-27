package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.UNSET;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ConsensusState
{
    public final boolean IsForward;
    public final String Chromosome;
    private final RefGenomeInterface mRefGenome;
    public byte[] Bases;
    public byte[] BaseQualities;
    public List<CigarElement> CigarElements;

    public int MinUnclippedPosStart;
    public int MaxUnclippedPosEnd;
    public int MinAlignedPosStart;
    public int MaxAlignedPosEnd;
    public int MapQuality;
    public int NumMutations;

    private ConsensusOutcome mOutcome;

    public ConsensusState(final boolean isForward, final String chromosome, final RefGenomeInterface refGenome)
    {
        IsForward = isForward;
        Chromosome = chromosome;
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

    public void addCigarElement(int length, final CigarOperator operator)
    {
        // combine with existing if a match on type
        if(IsForward)
        {
            int lastIndex = CigarElements.size() - 1;
            if(lastIndex >= 0 && CigarElements.get(lastIndex).getOperator() == operator)
                CigarElements.set(lastIndex, new CigarElement(CigarElements.get(lastIndex).getLength() + length, operator));
            else
                CigarElements.add(new CigarElement(length, operator));
        }
        else
        {
            int firstIndex = !CigarElements.isEmpty() ? 0 : -1;
            if(firstIndex >= 0 && CigarElements.get(firstIndex).getOperator() == operator)
                CigarElements.set(firstIndex, new CigarElement(CigarElements.get(firstIndex).getLength() + length, operator));
            else
                CigarElements.add(0, new CigarElement(length, operator));
        }
    }

    public void setNumMutations()
    {
        NumMutations = 0;
        String refBases = mRefGenome.getBaseString(Chromosome, MinAlignedPosStart, MaxAlignedPosEnd);
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
                    if(refBases.charAt(refBaseIndex + i) != (char) Bases[baseIndex + i])
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
