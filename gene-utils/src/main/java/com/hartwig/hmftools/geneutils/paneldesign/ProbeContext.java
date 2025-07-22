package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Encapsulates the constant fields so we don't have to pass them around everywhere when creating probe probes.
public record ProbeContext(
        TargetRegion targetRegion
)
{
    public Probe createProbe(final ChrBaseRegion probeRegion)
    {
        return new Probe(targetRegion, probeRegion);
    }

    public Probe createProbe(final BaseRegion probeRegion)
    {
        return createProbe(ChrBaseRegion.from(targetRegion.region().chromosome(), probeRegion));
    }
}
