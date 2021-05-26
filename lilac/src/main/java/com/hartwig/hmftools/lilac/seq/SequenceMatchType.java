package com.hartwig.hmftools.lilac.seq;

public enum SequenceMatchType
{
    FULL,
    PARTIAL,
    WILD,
    NONE;

    public boolean isBetter(final SequenceMatchType other)
    {
        return this.rank() < other.rank();
    }

    public int rank()
    {
        if(this == FULL) return 0;
        if(this == PARTIAL) return 1;
        if(this == WILD) return 2;
        return 3;
    }

}
