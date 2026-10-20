package com.hartwig.hmftools.viridian.profiling.support;

// For a virus genome: is it possibly in the sample, or is it noise?
public enum GenomeFilterStatus
{
    // Cleared both filters.
    PASS,
    // Not enough coverage.
    LOW_COVERAGE,
    // Passed coverage but not enough read votes.
    // Can occur if marginally similar to another genome, but that other genome gets much better alignments.
    LOW_VOTE_DENSITY
}
