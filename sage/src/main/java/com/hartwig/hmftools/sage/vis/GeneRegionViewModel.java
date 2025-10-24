package com.hartwig.hmftools.sage.vis;

import com.hartwig.hmftools.common.region.BaseRegion;

public interface GeneRegionViewModel
{
    record AminoAcidViewModel_(BaseRegion baseRegion, int aminoAcidPos, char ref, char alt) implements GeneRegionViewModel
    {
        public boolean matchesRef()
        {
            return ref == alt;
        }
    }

    record IntronicRegionViewModel(BaseRegion baseRegion) implements GeneRegionViewModel {}
}
