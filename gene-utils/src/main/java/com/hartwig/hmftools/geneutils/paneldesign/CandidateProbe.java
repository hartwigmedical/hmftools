package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public record CandidateProbe(
        // What caused the probe to be generated.
        ProbeSourceInfo source,
        // The region of the probe itself.
        ChrBaseRegion probeRegion,
        // The region we were originally interested to create a probe for. Just for output and debugging purposes.
        ChrBaseRegion targetRegion
)
{
    public CandidateProbe
    {
        if(!probeRegion.Chromosome.equals(targetRegion.Chromosome))
        {
            // Not expecting this to ever occur but should check it to be safe.
            throw new IllegalArgumentException("probeRegion and targetRegion should have same chromosome");
        }
    }
}
