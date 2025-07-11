package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Regions which was requested to be covered by probes but couldn't be due to probe selection criteria.
public record RejectedRegion(
        ChrBaseRegion baseRegion,
        ProbeSourceInfo source,
        String reason
)
{
    public static RejectedRegion fromTargetRegion(TargetRegion target, String reason)
    {
        return new RejectedRegion(target.region(), target.source(), reason);
    }
}
