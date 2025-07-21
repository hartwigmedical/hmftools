package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// A probe which is plausible for a particular target region, subject to future evaluation and filtering.
public record CandidateProbe(
        TargetRegion target,
        ChrBaseRegion probeRegion
)
{
    public CandidateProbe
    {
        if(!probeRegion.chromosome().equals(target.region().chromosome()))
        {
            // Not expecting this to ever occur but should check it to be safe.
            throw new IllegalArgumentException("Probe region and target region should have the same chromosome");
        }
    }
}
