package com.hartwig.hmftools.viridian.detection.support;

// Filter status of a virus genome, derived from read support.
public enum GenomeFilterStatus
{
    // Cleared both filters.
    CANDIDATE,
    // Not enough coverage.
    LOW_COVERAGE,
    // Passed coverage but not enough read votes.
    // Can occur if marginally similar to another genome, but that other genome gets much better alignments.
    LOW_VOTE_DENSITY
}
