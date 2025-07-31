package com.hartwig.hmftools.redux.bqr;

import org.checkerframework.checker.units.qual.N;

public class AltQualityCount
{
    public final byte Alt;
    public final byte Quality;
    public int PosStrandCount;
    public int NegStrandCount;

    public AltQualityCount(final byte alt, final byte quality, final boolean posStrand)
    {
        Alt = alt;
        Quality = quality;
        PosStrandCount = 0;
        NegStrandCount = 0;

        increment(posStrand);
    }

    public void increment(final boolean posStrand)
    {
        if(posStrand)
            ++PosStrandCount;
        else
            ++NegStrandCount;
    }

    public int totalCount() { return PosStrandCount + NegStrandCount; }

    public String toString()
    {
        return String.format("alt(%s) qual(%d) count(pos=%d neg=%d)",
                (char) Alt, (int) Quality, PosStrandCount, NegStrandCount);
    }
}
