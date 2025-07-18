package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Region which was requested to be (fully or partially) covered by probes.
// Target regions ignore overlaps. A target region may not be fully covered by probes if a subregion was already covered by other probes.
// Just for output and debugging purposes.
public record TargetRegion(
        // This region may or may not be meaningful depending on the probe source type and how the probes are selected.
        ChrBaseRegion region,
        TargetMetadata metadata
)
{
}
