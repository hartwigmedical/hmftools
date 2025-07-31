package com.hartwig.hmftools.panelbuilder;

import com.hartwig.hmftools.common.region.BasePosition;

// Amber heterozygous site data.
public record AmberSite(
        BasePosition position,
        double gnomadFreq
)
{
}
