package com.hartwig.hmftools.virusdetect;

import org.jetbrains.annotations.Nullable;

// Representative contig selection's outcome for one viral contig.
public record ContigSelectionResult(
        ContigSupport stats,
        ContigFilterStatus filterStatus,
        // Null if the contig was prefiltered and not considered as a candidate.
        @Nullable CandidateSelectionResult candidate
)
{
    public ViralContig contig()
    {
        return stats.contig();
    }
}
