package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Region which was requested to be covered by probes. Just for output and debugging purposes.
public record TargetRegion(
        ChrBaseRegion region,
        TargetMetadata metadata
)
{
}
