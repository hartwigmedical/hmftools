package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public record CustomRegion(
        ChrBaseRegion region,
        // Arbitrary descriptor for the user.
        String extraInfo
)
{
}
