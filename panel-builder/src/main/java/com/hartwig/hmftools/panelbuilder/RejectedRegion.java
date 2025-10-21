package com.hartwig.hmftools.panelbuilder;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// TODO: should this include rejected sequences for variant probes?
// Regions which was requested to be covered by probes but couldn't be due to probe selection criteria.
// Just for output and debugging purposes.
public record RejectedRegion(
        // The exact region that could not be covered.
        ChrBaseRegion region,
        TargetMetadata metadata,
        // TODO: the reason is not very useful. should just have the eval criteria
        // The reason the region could not be covered.
        String reason
)
{
    public RejectedRegion
    {
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
    }
}
