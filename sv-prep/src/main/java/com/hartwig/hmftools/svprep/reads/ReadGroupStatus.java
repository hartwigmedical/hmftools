package com.hartwig.hmftools.svprep.reads;

public enum ReadGroupStatus
{
    UNSET,
    INCOMPLETE, // missing a mate for a read (matching supplementary status)
    EXPECTED,
    PAIRED, // all non-supp reads have their mates
    SUPPLEMENTARY, // all reads are supplementaries
    COMPLETE, // has all reads
}
