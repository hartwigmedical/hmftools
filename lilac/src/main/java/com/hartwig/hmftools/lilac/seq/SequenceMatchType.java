package com.hartwig.hmftools.lilac.seq;

public enum SequenceMatchType
{
    FULL,
    WILD,
    NO_LOCI,
    MISMATCH;

    public boolean isBetter(final SequenceMatchType other)
    {
        return this.rank() < other.rank();
    }

    public int rank()
    {
        if(this == FULL) return 0;
        if(this == WILD) return 1;
        if(this == NO_LOCI) return 1;
        return 3;
    }

}
