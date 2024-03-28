package com.hartwig.hmftools.sage.bqr;

public class AltQualityCount
{
    public final byte Alt;
    public final byte Quality;
    public int Count;

    public AltQualityCount(final byte alt, final byte quality)
    {
        Alt = alt;
        Quality = quality;
        Count = 1;
    }

    public String toString()
    {
        return String.format("alt(%s) qual(%d) count(%d)", (char) Alt, (int) Quality, Count);
    }
}
