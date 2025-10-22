package com.hartwig.hmftools.sage.vis;

import com.hartwig.hmftools.common.region.BaseRegion;

public record AminoAcidViewModel(BaseRegion baseRegion, int aminoAcidPos, char ref, char alt)
{
    public boolean matchesRef()
    {
        return ref == alt;
    }
}
