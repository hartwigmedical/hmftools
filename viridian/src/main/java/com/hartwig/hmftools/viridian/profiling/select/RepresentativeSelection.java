package com.hartwig.hmftools.viridian.profiling.select;

import java.util.List;

import com.hartwig.hmftools.viridian.reference.VirusGenome;

// One main virus genome selected per oncology group.
public record RepresentativeSelection(
        VirusGenome representativeGenome,
        // Other genomes just for output purposes. Ordered by most votes first.
        List<VirusGenome> secondaryGenomes
)
{
}
