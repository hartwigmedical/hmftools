package com.hartwig.hmftools.virusdetect;

// Outcome of the viral contig selection prefilters.
public enum ContigFilterStatus
{
    // Cleared both prefilters.
    CANDIDATE,
    // Not enough coverage.
    LOW_COVERAGE,
    // Passed coverage but not enough read votes.
    // Can occur if marginally similar to another contig, but that other contig gets much better alignments.
    LOW_VOTE_DENSITY
}
