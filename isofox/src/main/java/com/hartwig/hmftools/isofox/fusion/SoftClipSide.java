package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.isofox.common.ReadRecord;

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

    public static SoftClipSide fromRead(final Cigar cigar)
    {
        int scLeft = cigar.isLeftClipped() ? cigar.getFirstCigarElement().getLength() : 0;
        int scRight = cigar.isRightClipped() ? cigar.getLastCigarElement().getLength() : 0;

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
