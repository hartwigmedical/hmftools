package com.hartwig.hmftools.common.variant;

import org.apache.logging.log4j.util.Strings;

public class MicrohomologyContext
{
    public final byte[] Bases;
    private final int LeftAlignedPosition;
    private final int Length;

    public MicrohomologyContext(final int leftAlignedPosition, final byte[] bases, final int length)
    {
        Bases = bases;
        LeftAlignedPosition = leftAlignedPosition;
        Length = length;
    }

    public int position()
    {
        return LeftAlignedPosition;
    }

    public byte[] readSequence()
    {
        return Bases;
    }

    public int length()
    {
        return Length;
    }

    @Override
    public String toString()
    {
        return Length > 0 ? new String(Bases, homologyIndex(), Length) : Strings.EMPTY;
    }

    public int homologyIndex()
    {
        return LeftAlignedPosition + 1;
    }
}
