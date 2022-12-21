package com.hartwig.hmftools.bamtools.markdups;

public enum FragmentStatus
{
    UNSET,
    NONE,
    UNCLEAR,
    PRIMARY,
    DUPLICATE,
    SUPPLEMENTARY;

    public boolean isResolved() { return this == PRIMARY || this == DUPLICATE || this == NONE; }
    public boolean isDuplicate() { return this == PRIMARY || this == DUPLICATE; }
}
