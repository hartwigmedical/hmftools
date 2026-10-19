package com.hartwig.hmftools.viridian.detection.select;

import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

public record RepresentativeGenome(
        VirusGenome genome,
        RepresentativeStatus status
)
{
    public OncologyGroup oncologyGroup()
    {
        return genome.oncologyGroup();
    }
}
