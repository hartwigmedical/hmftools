package com.hartwig.hmftools.virusdetect;

// Coarse outcome of representative contig selection for an oncology group.
public enum OncologyGroupResolution
{
    // Nothing passed the prefilter, so no representative was sought.
    NO_CANDIDATES,
    // Representative found.
    RESOLVED,
    // No representative possible.
    UNRESOLVED
}
