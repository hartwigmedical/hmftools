package com.hartwig.hmftools.esvee.assembly.types;

public enum AssemblyOutcome
{
    UNSET,
    NO_LINK,
    LINKED, // 2 assemblies linked in a standard SV
    SECONDARY, // an assembly linked to an assembly which was primarily linked to another
    DUP_BRANCHED, // an
    LOCAL_INDEL, // an assembly linked to a local ref-genome sequence as a DEL or DUP
    REMOTE_REGION, // an assembly matched one or more of its remote region read sequences
    SUPP_ONLY; // the assmebly comprised supplementary reads only

    public boolean isDuplicate() { return this == DUP_BRANCHED; }
}
