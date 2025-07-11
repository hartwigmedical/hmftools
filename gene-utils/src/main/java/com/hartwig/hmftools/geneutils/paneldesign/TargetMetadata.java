package com.hartwig.hmftools.geneutils.paneldesign;

// Metadata information about a target region that we want to cover with probes.
public record TargetMetadata(
        TargetRegionType type,
        String extra
)
{
}
