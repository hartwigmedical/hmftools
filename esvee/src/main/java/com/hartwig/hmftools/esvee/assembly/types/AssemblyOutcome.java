package com.hartwig.hmftools.esvee.assembly.types;

public enum AssemblyOutcome
{
    LOCAL_INDEL, // an assembly linked to a local ref-genome sequence as a DEL or DUP
    LINKED, // 2 assemblies linked in a standard SV
    REMOTE_LINK, // an assembly matched one or more of its remote region read sequences
    REMOTE_REGION, // the remote region read matched from a standard assembly extension
    DUP_BRANCHED, // an assembly is branched from ref-base differences to allow it to make multiple links
    SECONDARY, // an assembly linked to an assembly which was primarily linked to another
    SUPP_ONLY, // the assmebly comprised supplementary reads only
    NO_LINK,
    DECOY, // matched a decoy
    UNSET;
}
