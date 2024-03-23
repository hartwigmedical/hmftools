package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipped;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static htsjdk.samtools.CigarOperator.S;

import htsjdk.samtools.Cigar;

public class ClippedSide
{
    public final int Side;
    public final int Length;
    public final boolean IsSoft;

    public ClippedSide(final int side, final int length, final boolean isSoft)
    {
        Side = side;
        Length = length;
        IsSoft = isSoft;
    }

    public static ClippedSide fromCigar(final Cigar cigar, boolean allowHardClips)
    {
        return from(
                cigar,
                allowHardClips ? cigar.isLeftClipped() : leftSoftClipped(cigar),
                allowHardClips ? cigar.isRightClipped() : rightSoftClipped(cigar));
    }

    public static ClippedSide from(final Cigar cigar, boolean isLeftClipped, boolean isRightClipped)
    {
        int scLeft = isLeftClipped ? cigar.getFirstCigarElement().getLength() : 0;
        int scRight = isRightClipped ? cigar.getLastCigarElement().getLength() : 0;

        if(scLeft == 0 && scRight == 0)
            return null;

        if(scLeft > 0 && scRight > 0)
        {
            if(scLeft >= scRight)
                return new ClippedSide(SE_START, scLeft, cigar.getFirstCigarElement().getOperator() == S);
            else
                return new ClippedSide(SE_END, scRight, cigar.getLastCigarElement().getOperator() == S);
        }
        else if(scLeft > 0)
        {
            return new ClippedSide(SE_START, scLeft, cigar.getFirstCigarElement().getOperator() == S);
        }
        else
        {
            return new ClippedSide(SE_END, scRight, cigar.getLastCigarElement().getOperator() == S);
        }
    }

    public boolean isLeft() { return Side == SE_START; }
    public boolean isRight() { return Side == SE_END; }
}
