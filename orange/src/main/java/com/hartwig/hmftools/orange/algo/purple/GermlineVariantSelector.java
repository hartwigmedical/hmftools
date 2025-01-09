package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.Nullable;

final class GermlineVariantSelector
{
    @Nullable
    public static List<PurpleVariant> selectInterestingUnreportedVariants(@Nullable List<PurpleVariant> allGermlineVariants)
    {
        if(allGermlineVariants == null)
        {
            return null;
        }

        List<PurpleVariant> filtered = Lists.newArrayList();
        for(PurpleVariant variant : allGermlineVariants)
        {
            if(!variant.reported())
            {
                boolean isHotspot = variant.hotspot() == HotspotType.HOTSPOT;

                // TODO: Add pathogenic variants that were not reported
                // TODO: Add variants with conflicting evidence in ClinVar
                if(isHotspot)
                {
                    filtered.add(variant);
                }
            }
        }
        return filtered;
    }
}
