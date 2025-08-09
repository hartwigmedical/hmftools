package com.hartwig.hmftools.sage.seqtech;

public class RefMask
{
    public final int PosStart;
    public final int PosEnd;
    public final byte BaseMask;

    public RefMask(int posStart, int posEnd, byte baseMask)
    {
        PosStart = posStart;
        PosEnd = posEnd;
        BaseMask = baseMask;
    }

    // TODO: Remove
    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
            return true;

        if(!(o instanceof RefMask))
            return false;

        final RefMask refMask = (RefMask) o;
        return PosStart == refMask.PosStart && PosEnd == refMask.PosEnd && BaseMask == refMask.BaseMask;
    }

    @Override
    public int hashCode()
    {
        int hash = PosStart;
        hash = PosEnd + 31 * hash;
        hash = BaseMask + 31 * hash;
        return hash;
    }
}
