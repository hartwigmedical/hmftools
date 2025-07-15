package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Encapsulates the constant fields so we don't have to pass them around everywhere when creating candidate probes.
public record ProbeFactory(
        TargetRegion targetRegion
)
{
    public CandidateProbe create(final ChrBaseRegion probeRegion)
    {
        return new CandidateProbe(targetRegion, probeRegion);
    }
}
