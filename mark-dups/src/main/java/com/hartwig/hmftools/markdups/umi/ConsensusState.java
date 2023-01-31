package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_CONSENSUS_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.UNSET;
import static com.hartwig.hmftools.markdups.umi.UmiConfig.READ_ID_DELIM;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ConsensusState
{
    public final boolean IsForward;
    public final String Chromosome;
    public byte[] Bases;
    public byte[] BaseQualities;
    public List<CigarElement> CigarElements;

    public int MinUnclippedPosStart;
    public int MaxUnclippedPosEnd;
    public int MinAlignedPosStart;
    public int MaxAlignedPosEnd;
    public int MapQuality;

    private ConsensusOutcome mOutcome;

    public ConsensusState(final boolean isForward, final String chromosome)
    {
        IsForward = isForward;
        Chromosome = chromosome;
        Bases = null;
        BaseQualities = null;
        CigarElements = Lists.newArrayList();

        MinUnclippedPosStart = 0;
        MaxUnclippedPosEnd = 0;
        MinAlignedPosStart = 0;
        MaxAlignedPosEnd = 0;
        MapQuality = 0;

        mOutcome = UNSET;
    }

    public ConsensusOutcome outcome() { return mOutcome; }
    public void setOutcome(final ConsensusOutcome outcome) { mOutcome = outcome; }

    void setBaseLength(int baseLength)
    {
        Bases = new byte[baseLength];
        BaseQualities = new byte[baseLength];
    }

    void expandBaseLength(int baseLength)
    {
        if(Bases == null)
        {
            setBaseLength(baseLength);
            return;
        }

        if(IsForward)
        {
            Bases = Arrays.copyOf(Bases, Bases.length + baseLength);
            BaseQualities = Arrays.copyOf(BaseQualities, BaseQualities.length + baseLength);
        }
        else
        {
            byte[] bases = new byte[Bases.length + baseLength];
            byte[] baseQuals = new byte[Bases.length + baseLength];

            for(int i = 0; i < Bases.length; ++i)
            {
                bases[i + baseLength] = Bases[i];
                baseQuals[i + baseLength] = BaseQualities[i];
            }

            Bases = bases;
            BaseQualities = baseQuals;
        }
    }

    void setBoundaries(final SAMRecord read)
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
}
