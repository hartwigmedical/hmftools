package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public record ProbeCandidate(
        // What caused the probe to be generated.
        ProbeSourceInfo source,
        // The region of the probe itself.
        ChrBaseRegion probeRegion,
        // TODO: needed?
        // The region we were originally interested to create a probe for.
        // In the best case, it's contained within probeRegion. Depending on probe selection criteria, it may only mostly overlap.
        ChrBaseRegion targetRegion
)
{
}
