package com.hartwig.hmftools.common.samtools;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;

public class SoftClipSide
{
    public final int Side;
    public final int Length;

    public SoftClipSide(final int side, final int length)
    {
        Side = side;
        Length = length;
    }

    public static SoftClipSide fromCigar(final Cigar cigar) { return from(cigar, cigar.isLeftClipped(), cigar.isRightClipped()); }

    public static SoftClipSide from(final Cigar cigar, boolean isLeftClipped, boolean isRightClipped)
    {
        int scLeft = isLeftClipped ? cigar.getFirstCigarElement().getLength() : 0;
        int scRight = isRightClipped ? cigar.getLastCigarElement().getLength() : 0;

        if(scLeft == 0 && scRight == 0)
            return null;

        if(scLeft > 0 && scRight > 0)
        {
            if(scLeft >= scRight)
                return new SoftClipSide(SE_START, scLeft);
            else
                return new SoftClipSide(SE_END, scRight);
        }
        else if(scLeft > 0)
        {
            return new SoftClipSide(SE_START, scLeft);
        }
        else
        {
            return new SoftClipSide(SE_END, scRight);
        }
    }

    public boolean isLeft() { return Side == SE_START; }
    public boolean isRight() { return Side == SE_END; }
}
