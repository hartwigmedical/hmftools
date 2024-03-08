package com.hartwig.hmftools.esvee.common;

public enum AssemblyOutcome
{
    UNSET,
    NO_LINK,
    SECONDARY,
    DUP_SPLIT,
    DUP_BRANCHED,
    SUPP_ONLY,
    LINKED;

    public boolean isDuplicate() { return this == DUP_BRANCHED || this == DUP_SPLIT; }
}
