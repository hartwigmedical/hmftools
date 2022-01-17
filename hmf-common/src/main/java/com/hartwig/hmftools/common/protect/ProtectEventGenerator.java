package com.hartwig.hmftools.common.protect;

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

        return variant.canonicalEffect();
    }

    @NotNull
    public static String fusionEvent(@NotNull LinxFusion fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd() + " fusion";
    }

}
