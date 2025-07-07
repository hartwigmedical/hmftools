package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.region.BasePosition;

// Amber heterozygous site data.
public record AmberSite(
        BasePosition position,
        double gnomadFreq,
        double mappability,
        double gcRatio
)
{
}
