package com.hartwig.hmftools.common.protect;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

public final class ProtectEventGenerator {

    static final String UPSTREAM_GENE_VARIANT = "upstream_gene_variant";

    private ProtectEventGenerator() {
    }

    @NotNull
    public static String variantEvent(@NotNull Variant variant) {
        return variant(variant.canonicalHgvsProteinImpact(), variant.canonicalHgvsCodingImpact(), variant.canonicalEffect());
    }

    @NotNull
    public static String variantEventNonCanonical(@NotNull String otherReportedEffects) {
        String[] parts = otherReportedEffects.split("\\|");
        return variant(parts[2], parts[1], parts[3]);
    }

    @NotNull
    private static String variant(@NotNull String protein, @NotNull String coding, @NotNull String effect) {
        if (!protein.isEmpty()) {
            return protein;
        }

        if (!coding.isEmpty()) {
            return coding;
        }

        if (effect.equals(UPSTREAM_GENE_VARIANT)) {
            return "upstream";
        }

        return effect;
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
