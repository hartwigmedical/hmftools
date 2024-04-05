package com.hartwig.hmftools.esvee.assembly.types;

public enum AssemblyOutcome
{
    UNSET,
    NO_LINK,
    LINKED,
    SHORT_INDEL,
    SECONDARY,
    DUP_SPLIT,
    DUP_BRANCHED,
    REMOTE_REF,
    SUPP_ONLY;

    public boolean isDuplicate() { return this == DUP_BRANCHED || this == DUP_SPLIT; }
}
