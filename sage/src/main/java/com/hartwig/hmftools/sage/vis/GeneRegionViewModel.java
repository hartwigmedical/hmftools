package com.hartwig.hmftools.sage.vis;

import com.hartwig.hmftools.common.region.BaseRegion;

public interface GeneRegionViewModel
{
    record AminoAcidViewModel(BaseRegion baseRegion, int aminoAcidPos, char ref, char alt, int insertAfterLength) implements GeneRegionViewModel
    {
        public boolean matchesRef()
        {
            return ref == alt;
        }
    }

    record NonCodingExonicRegionViewModel(BaseRegion baseRegion) implements GeneRegionViewModel {}
    record IntronicRegionViewModel(BaseRegion baseRegion) implements GeneRegionViewModel {}
}
