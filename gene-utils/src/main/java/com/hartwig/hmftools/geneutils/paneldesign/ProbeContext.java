package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Encapsulates the constant fields so we don't have to pass them around everywhere when creating probe probes.
public record ProbeContext(
        TargetMetadata metadata
)
{
    public Probe createProbe(final ChrBaseRegion region)
    {
        return new Probe(region, metadata);
    }
}
