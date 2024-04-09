package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

public class IndelCoords
{
    public final int PosStart;
    public final int PosEnd;
    public final int Length;

    private String mInsertedBases;

    public IndelCoords(final int posStart, final int posEnd, final int length)
    {
        PosStart = posStart;
        PosEnd = posEnd;
        Length = length;
        mInsertedBases = null;
    }

    public boolean isInsert() { return PosEnd == PosStart + 1; }
    public boolean isDelete() { return !isInsert(); }

    public String insertedBases() { return mInsertedBases != null ? mInsertedBases : ""; }
    public void setInsertedBases(final String bases) { mInsertedBases = bases; }

    public boolean matchesJunction(int position, byte orientation)
    {
        return orientation == POS_ORIENT ? PosStart == position : PosEnd == position;
    }

    public boolean matches(final IndelCoords other)
    {
        return PosStart == other.PosStart && PosEnd == other.PosEnd && Length == other.Length;
    }

    public String toString()
    {
        return format("%s(%d - %d) len(%d)", isDelete() ? "delete" : "insert", PosStart, PosEnd, Length);
    }
}
