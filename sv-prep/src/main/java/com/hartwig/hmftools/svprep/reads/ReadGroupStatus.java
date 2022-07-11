package com.hartwig.hmftools.svprep.reads;

public enum ReadGroupStatus
{
    UNSET,
    INCOMPLETE,
    PARTIAL, // has all reads for initial partial, others are remote
    COMPLETE;
}
