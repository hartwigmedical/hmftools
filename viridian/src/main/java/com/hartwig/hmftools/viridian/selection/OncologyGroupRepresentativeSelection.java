package com.hartwig.hmftools.viridian.selection;

import java.util.List;

import com.hartwig.hmftools.viridian.detection.contig_support.ContigSupport;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.jetbrains.annotations.Nullable;

// Representative contig selection's outcome for one oncology group.
// Includes all the contigs which had any alignments.
public record OncologyGroupRepresentativeSelection(
        OncologyGroup oncologyGroup,
        OncologyGroupOutcome outcome,
        List<RepresentativeContigCandidate> candidates,
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
                .map(RepresentativeContigCandidate::contig)
                .findFirst()
                .orElse(null);
    }
}
