package com.hartwig.hmftools.esvee.assembly;

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
