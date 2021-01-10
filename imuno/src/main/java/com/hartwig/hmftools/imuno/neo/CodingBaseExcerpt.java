package com.hartwig.hmftools.imuno.neo;

public class CodingBaseExcerpt
{
    public final String Bases;
    public final int[] Positions;

    public CodingBaseExcerpt(final String bases, final int upPos, final int downPos)
    {
        Bases = bases;
        Positions = new int[] { upPos, downPos};
    }
}
