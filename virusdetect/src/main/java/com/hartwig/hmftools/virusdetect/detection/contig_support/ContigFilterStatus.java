package com.hartwig.hmftools.virusdetect.detection.contig_support;

// Filter status of a viral genome contig, derived from read support.
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
