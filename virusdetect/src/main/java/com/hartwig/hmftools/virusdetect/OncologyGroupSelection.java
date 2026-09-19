package com.hartwig.hmftools.virusdetect;

import java.util.List;

import org.jetbrains.annotations.Nullable;

// Representative contig selection's outcome for one oncology group.
// Includes all the contigs which had any alignments.
public record OncologyGroupSelection(
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

    public int votesRank(ViralContig contig)
    {
        for(int i = 0; i < candidates.size(); ++i)
        {
            if(candidates.get(i).contig().equals(contig))
            {
                return i + 1;
            }
        }
        throw new IllegalArgumentException("Contig was not a selection candidate: " + contig.name());
    }
}
