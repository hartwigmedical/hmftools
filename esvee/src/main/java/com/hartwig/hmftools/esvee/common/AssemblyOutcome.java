package com.hartwig.hmftools.esvee.common;

public enum AssemblyOutcome
{
    UNSET,
    NO_LINK,
    LINKED,
    SHORT_INDEL,
    SECONDARY,
    DUP_SPLIT,
    DUP_BRANCHED,
    SUPP_ONLY;

    public boolean isDuplicate() { return this == DUP_BRANCHED || this == DUP_SPLIT; }
}
