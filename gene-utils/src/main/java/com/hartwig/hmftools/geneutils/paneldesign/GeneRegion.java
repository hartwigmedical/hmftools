package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.BaseRegion;

public record GeneRegion(
        GeneInfo gene,
        GeneRegionType type,
        BaseRegion region
)
{
}
