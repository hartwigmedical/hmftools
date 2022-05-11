package com.hartwig.hmftools.common.protect;

import com.hartwig.hmftools.common.protect.variant.OtherEffectsInterpreter;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

public final class ProtectEventGenerator {

    static final String UPSTREAM_GENE_VARIANT = "upstream_gene_variant";

    private ProtectEventGenerator() {
    }

    @NotNull
    public static String variantEvent(@NotNull Variant variant) {
        if (variant instanceof ReportableVariant) {
            return reportableVariantEvent((ReportableVariant) variant);
        } else {
            return canonicalVariantEvent(variant);
        }
    }

    @NotNull
    private static String reportableVariantEvent(@NotNull ReportableVariant reportableVariant) {
        return reportableVariant.isCanonical() ? canonicalVariantEvent(reportableVariant) : nonCanonicalVariantEvent(reportableVariant);
    }

    @NotNull
    private static String canonicalVariantEvent(@NotNull Variant variant) {
        return toVariantEvent(variant.canonicalHgvsProteinImpact(), variant.canonicalHgvsCodingImpact(), variant.canonicalEffect());
    }

    @NotNull
    private static String nonCanonicalVariantEvent(@NotNull ReportableVariant variant) {
        return toVariantEvent(OtherEffectsInterpreter.hgvsProteinImpact(variant.otherReportedEffects()),
                OtherEffectsInterpreter.hgvsCodingImpact(variant.otherReportedEffects()),
                OtherEffectsInterpreter.effect(variant.otherReportedEffects()));
    }

    @NotNull
    private static String toVariantEvent(@NotNull String protein, @NotNull String coding, @NotNull String effect) {
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
