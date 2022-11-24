package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class GermlineVariantSelector {

    private GermlineVariantSelector() {
    }

    @Nullable
    public static List<ReportableVariant> selectInterestingUnreportedVariants(@Nullable List<SomaticVariant> allGermlineVariants) {
        if (allGermlineVariants == null) {
            return null;
        }

        List<ReportableVariant> filtered = Lists.newArrayList();
        for (SomaticVariant variant : allGermlineVariants) {
            if (!variant.reported()) {
                boolean isHotspot = variant.isHotspot();

                // TODO: Add pathogenic variants that were not reported
                // TODO: Add variants with conflicting evidence
                if (isHotspot) {
                    filtered.add(toReportable(variant));
                }
            }
        }
        return filtered;
    }

    @NotNull
    private static ReportableVariant toReportable(@NotNull SomaticVariant variant) {
        return ReportableVariantFactory.fromVariant(variant, ReportableVariantSource.GERMLINE)
                .driverLikelihood(Double.NaN)
                .transcript(variant.canonicalTranscript())
                .isCanonical(true)
                .build();
    }
}
