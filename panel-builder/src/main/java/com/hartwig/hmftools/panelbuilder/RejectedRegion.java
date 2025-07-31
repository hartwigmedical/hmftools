package com.hartwig.hmftools.panelbuilder;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Regions which was requested to be covered by probes but couldn't be due to probe selection criteria.
public record RejectedRegion(
        // The exact region that could not be covered.
        ChrBaseRegion region,
        TargetMetadata metadata,
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

    // Convenience method for creating a rejected region from an entire target region.
    public static RejectedRegion rejectTarget(TargetRegion target, String reason)
    {
        return new RejectedRegion(target.region(), target.metadata(), reason);
    }
}
