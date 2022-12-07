package com.hartwig.hmftools.sage.evidence;

public enum SyncFragmentType
{
    COMBINED,
    NO_OVERLAP,
    CIGAR_MISMATCH,
    NO_OVERLAP_CIGAR_DIFF,
    BASE_MISMATCH,
    EXCEPTION;

    public boolean processSeparately()
    {
        return this == NO_OVERLAP || this == NO_OVERLAP_CIGAR_DIFF;
    }
}
