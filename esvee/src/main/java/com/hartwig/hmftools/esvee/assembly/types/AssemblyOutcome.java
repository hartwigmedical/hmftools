package com.hartwig.hmftools.esvee.assembly.types;

public enum AssemblyOutcome
{
    LOCAL_INDEL, // an assembly linked to a local ref-genome sequence as a DEL or DUP
    LINKED, // 2 assemblies linked in a standard SV
    DUP_BRANCHED, // an assembly is branched from ref-base differences to allow it to make multiple links
    SECONDARY, // an assembly linked to an assembly which was primarily linked to another
    SUPP_ONLY, // the assembly comprised supplementary reads only
    NO_LINK,
    DECOY, // matched a decoy
    UNSET;
}
