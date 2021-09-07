package com.hartwig.hmftools.lilac.fragment;

public enum FragmentScope
{
    BASE_QUAL_FILTERED,
    HLA_Y,
    NO_HET_LOCI,
    UNMATCHED_AMINO_ACID,
    WILD_ONLY,
    CANDIDATE,
    SOLUTION,
    UNSET;

    public boolean isUnmatched()
    {
        return this == CANDIDATE || this == UNMATCHED_AMINO_ACID || this == WILD_ONLY;
    }
}
