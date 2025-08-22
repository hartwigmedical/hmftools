package com.hartwig.hmftools.panelbuilder;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Region which was requested to be (fully or partially) covered by probes.
// Target regions ignore overlaps. A target region may not be fully covered by probes if a subregion was already covered by other probes.
// Just for output and debugging purposes.
public record TargetRegion(
        ChrBaseRegion region,
        TargetMetadata metadata
)
{
    public TargetRegion
    {
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
    }
}
