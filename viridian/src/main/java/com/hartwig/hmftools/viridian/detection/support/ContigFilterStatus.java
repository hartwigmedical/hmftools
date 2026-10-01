package com.hartwig.hmftools.viridian.detection.support;

// Filter status of a virus genome contig, derived from read support.
public enum ContigFilterStatus
{
    // Cleared both filters.
    CANDIDATE,
    // Not enough coverage.
    LOW_COVERAGE,
    // Passed coverage but not enough read votes.
    // Can occur if marginally similar to another contig, but that other contig gets much better alignments.
    LOW_VOTE_DENSITY
}
