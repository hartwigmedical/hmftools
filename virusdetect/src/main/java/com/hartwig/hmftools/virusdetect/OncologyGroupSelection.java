package com.hartwig.hmftools.virusdetect;

import java.util.List;

import org.jetbrains.annotations.Nullable;

// Representative contig selection's outcome for one oncology group.
// Includes all the contigs which had any alignments.
public record OncologyGroupSelection(
        OncologyGroup oncologyGroup,
        OncologyGroupOutcome outcome,
        List<ContigSelectionResult> contigs
)
{
    public OncologyGroupResolution resolution()
    {
        return outcome.resolution();
    }

    @Nullable
    public ViralContig representative()
    {
        return contigs.stream()
                .filter(contig -> contig.candidate() != null && contig.candidate().role() == ContigRole.REPRESENTATIVE)
                .map(ContigSelectionResult::contig)
                .findFirst()
                .orElse(null);
    }

    public List<ContigSelectionResult> candidates()
    {
        return contigs.stream().filter(contig -> contig.candidate() != null).toList();
    }
}
