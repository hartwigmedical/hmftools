package com.hartwig.hmftools.peach;

public enum Zygosity
{
    HOMOZYGOUS,
    HETEROZYGOUS;

    public static Zygosity fromInt(int count)
    {
        if (count == 2)
            return Zygosity.HOMOZYGOUS;
        else if (count == 1)
            return Zygosity.HETEROZYGOUS;
        else
            throw new RuntimeException(String.format("Cannot construct zygosity from count %d", count));
    }

    public int getVariantCount()
    {
        if (this == HOMOZYGOUS)
            return 2;
        else if (this == HETEROZYGOUS)
            return 1;
        else
            throw new RuntimeException(String.format("Cannot get count from zygosity %s", this));
    }
}
