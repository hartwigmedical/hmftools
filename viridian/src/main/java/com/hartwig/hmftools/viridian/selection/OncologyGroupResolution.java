package com.hartwig.hmftools.viridian.selection;

// Coarse outcome of representative contig selection for an oncology group.
public enum OncologyGroupResolution
{
    // Nothing passed the support filter, so no representative was sought.
    NO_CANDIDATES,
    // Representative found.
    RESOLVED,
    // No representative possible.
    UNRESOLVED
}
