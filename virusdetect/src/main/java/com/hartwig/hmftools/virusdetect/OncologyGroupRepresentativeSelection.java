package com.hartwig.hmftools.virusdetect;

import java.util.List;

import org.jetbrains.annotations.Nullable;

// Representative contig selection's outcome for one oncology group.
// Includes all the contigs which had any alignments.
public record OncologyGroupRepresentativeSelection(
        OncologyGroup oncologyGroup,
        OncologyGroupOutcome outcome,
        List<RepresentativeCandidate> candidates,
        // List of contigs which were filtered so not candidates for representative.
        List<ContigSupport> rejected
)
{
    public OncologyGroupResolution resolution()
    {
        return outcome.resolution();
    }

    @Nullable
    public ViralContig representative()
    {
        return candidates.stream()
                .filter(candidate -> candidate.role() == ContigRole.REPRESENTATIVE)
                .map(RepresentativeCandidate::contig)
                .findFirst()
                .orElse(null);
    }
}
