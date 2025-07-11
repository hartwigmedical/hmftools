package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public record CandidateProbe(
        TargetRegion target,
        ChrBaseRegion probeRegion
)
{
    public CandidateProbe
    {
        if(!probeRegion.Chromosome.equals(target.region().Chromosome))
        {
            // Not expecting this to ever occur but should check it to be safe.
            throw new IllegalArgumentException("Probe region and target region should have the same chromosome");
        }
    }
}
