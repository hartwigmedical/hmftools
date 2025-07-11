package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Regions which was requested to be covered by probes but couldn't be due to probe selection criteria.
public record RejectedRegion(
        // The exact region that could not be covered.
        ChrBaseRegion baseRegion,
        // The full region we were interested in covering.
        TargetRegion target,
        // The reason the region could not be covered.
        String reason
)
{
    public static RejectedRegion fromTargetRegion(TargetRegion target, String reason)
    {
        return new RejectedRegion(target.region(), target, reason);
    }
}
