package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class GermlineVariantSelector {

    private GermlineVariantSelector() {
    }

    @NotNull
    public static List<ReportableVariant> selectNonDrivers(@NotNull List<SomaticVariant> unreportedVariants) {
        List<ReportableVariant> filtered = Lists.newArrayList();
        for (SomaticVariant variant : unreportedVariants) {
            boolean isHotspot = variant.isHotspot();

            // TODO: Add pathogenic variants that were not reported
            // TODO: Add variants with conflicting evidence
            if (isHotspot) {
                filtered.add(toReportable(variant));
            }
        } return filtered;
    }

    @NotNull
    private static ReportableVariant toReportable(@NotNull SomaticVariant variant) {
        return ReportableVariantFactory.fromVariant(variant, ReportableVariantSource.GERMLINE).driverLikelihood(Double.NaN).build();
    }
}
