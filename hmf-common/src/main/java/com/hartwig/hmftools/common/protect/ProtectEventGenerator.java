package com.hartwig.hmftools.common.protect;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

public final class ProtectEventGenerator {

    private ProtectEventGenerator() {
    }

    @NotNull
    public static String variantEvent(@NotNull Variant variant) {
        String protein = variant.canonicalHgvsProteinImpact();
        if (!protein.isEmpty()) {
            return protein;
        }

        String coding = variant.canonicalHgvsCodingImpact();
        if (!coding.isEmpty()) {
            return coding;
        }

        // TODO: Also format "upstream" here, similar to VariantUtil.
        return variant.canonicalEffect();
    }

    @NotNull
    public static String variantEventNonCanonical(@NotNull String otherReportedEffects) {
        String protein = otherReportedEffects.split("\\|")[2];
        if (!protein.isEmpty()) {
            return protein;
        }

        String coding = otherReportedEffects.split("\\|")[1];
        if (!coding.isEmpty()) {
            return coding;
        }

        return otherReportedEffects.split("\\|")[3];
    }

    @NotNull
    public static String copyNumberEvent(@NotNull ReportableGainLoss gainLoss) {
        return gainLoss.interpretation().display();
    }

    @NotNull
    public static String fusionEvent(@NotNull LinxFusion fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd() + " fusion";
    }

}
