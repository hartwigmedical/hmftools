package com.hartwig.hmftools.redux.consensus;

public enum ConsensusOutcome
{
    UNSET,
    ALIGNMENT_ONLY,
    INDEL_MATCH,
    INDEL_MISMATCH,
    INDEL_FAIL,
    SUPPLEMENTARY;

    public boolean valid() { return this != INDEL_FAIL; }
}
