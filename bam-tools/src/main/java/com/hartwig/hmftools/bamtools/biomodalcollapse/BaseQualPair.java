package com.hartwig.hmftools.bamtools.biomodalcollapse;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

public class BaseQualPair
{
    public final byte Base;
    public final int Qual;

    public BaseQualPair(byte base, int qual)
    {
        Base = base;
        Qual = qual;
    }

    public BaseQualPair complementBase()
    {
        return new BaseQualPair(swapDnaBase(Base), Qual);
    }

    @Override
    public String toString()
    {
        return format("base(%s) qual(%d)", (char) Base, Qual);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof BaseQualPair))
        {
            return false;
        }
        final BaseQualPair that = (BaseQualPair) o;
        return Base == that.Base && Qual == that.Qual;
    }

    @Override
    public int hashCode()
    {
        return (int) Base + 31 * Qual;
    }
}
