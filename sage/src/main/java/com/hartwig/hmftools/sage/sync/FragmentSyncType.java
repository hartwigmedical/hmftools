package com.hartwig.hmftools.sage.sync;

public enum FragmentSyncType
{
    COMBINED,
    INVERSION,
    NO_OVERLAP,
    CIGAR_MISMATCH,
    BASE_MISMATCH,
    EXCEPTION;

    public boolean processSeparately()
    {
        return this == NO_OVERLAP || this == INVERSION;
    }
}
