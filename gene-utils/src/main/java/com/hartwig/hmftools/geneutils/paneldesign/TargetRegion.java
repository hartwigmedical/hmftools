package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Region which was requested to be (fully or partially) covered by probes. Just for output and debugging purposes.
public record TargetRegion(
        // This region may or may not be meaningful depending on the probe source type and how the probes are selected.
        ChrBaseRegion region,
        TargetMetadata metadata
)
{
}
