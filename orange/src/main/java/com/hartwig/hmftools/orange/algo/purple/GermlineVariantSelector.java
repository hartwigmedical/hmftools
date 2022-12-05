package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.variant.Hotspot;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.Nullable;

final class GermlineVariantSelector {

    private GermlineVariantSelector() {
    }

    @Nullable
    public static List<PurpleVariant> selectInterestingUnreportedVariants(@Nullable List<PurpleVariant> allGermlineVariants) {
        if (allGermlineVariants == null) {
            return null;
        }

        List<PurpleVariant> filtered = Lists.newArrayList();
        for (PurpleVariant variant : allGermlineVariants) {
            if (!variant.reported()) {
                boolean isHotspot = variant.hotspot() == Hotspot.HOTSPOT;

                // TODO: Add pathogenic variants that were not reported
                // TODO: Add variants with conflicting evidence
                if (isHotspot) {
                    filtered.add(variant);
                }
            }
        }
        return filtered;
    }
}
