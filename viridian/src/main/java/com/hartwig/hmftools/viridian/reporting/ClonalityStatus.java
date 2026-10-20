package com.hartwig.hmftools.viridian.reporting;

// The status of the virus clonality check for reporting, which can fail due to Purple failure.
public enum ClonalityStatus
{
    CLONAL,
    NOT_CLONAL,
    // Not present in sample reads.
    NOT_PROFILED,
    // Purple data was unusable.
    PURPLE_FAIL_NO_TUMOR,
    PURPLE_FAIL_CONTAMINATION
}
