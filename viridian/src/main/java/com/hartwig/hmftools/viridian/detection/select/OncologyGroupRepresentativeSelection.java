package com.hartwig.hmftools.viridian.detection.select;

import java.util.List;

import com.hartwig.hmftools.viridian.detection.support.GenomeSupport;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.jetbrains.annotations.Nullable;

// Representative selection's outcome for one oncology group.
// Includes all the virus genomes which had any alignments.
public record OncologyGroupRepresentativeSelection(
        OncologyGroup oncologyGroup,
        OncologyGroupOutcome outcome,
        List<RepresentativeGenomeCandidate> candidates,
        // List of virus genomes which were filtered so not candidates for representative.
        List<GenomeSupport> rejected
)
{
    public OncologyGroupResolution resolution()
    {
        return outcome.resolution();
    }

    @Nullable
    public VirusGenome representative()
    {
        return candidates.stream()
                .filter(candidate -> candidate.role() == GenomeRole.REPRESENTATIVE)
                .map(RepresentativeGenomeCandidate::genome)
                .findFirst()
                .orElse(null);
    }
}
