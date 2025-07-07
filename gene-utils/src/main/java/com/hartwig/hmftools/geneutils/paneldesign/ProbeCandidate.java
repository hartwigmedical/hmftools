package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public record ProbeCandidate(
        // What caused the probe to be generated.
        ProbeSourceInfo source,
        // The region of the probe itself.
        ChrBaseRegion probeRegion,
        // The region we were originally interested to create a probe for.
        // Ideally, contained within probeRegion. In less ideal cases it might partially overlap.
        ChrBaseRegion targetRegion
)
{
}
