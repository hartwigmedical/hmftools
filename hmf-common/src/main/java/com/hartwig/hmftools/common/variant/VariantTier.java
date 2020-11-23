package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public enum VariantTier {
    HOTSPOT,
    PANEL,
    HIGH_CONFIDENCE,
    LOW_CONFIDENCE,
    UNKNOWN;

    public static final String TIER = "TIER";

    @NotNull
    public static VariantTier fromContext(@NotNull final VariantContext context) {
        return fromString(context.getAttributeAsString(TIER, UNKNOWN.toString()));
    }

    @NotNull
    public static VariantTier fromString(@NotNull final String string) {
        switch (string) {
            case "HOTSPOT":
                return HOTSPOT;
            case "PANEL":
                return PANEL;
            case "HIGH_CONFIDENCE":
                return HIGH_CONFIDENCE;
            case "LOW_CONFIDENCE":
                return LOW_CONFIDENCE;
            default:
                return UNKNOWN;
        }
    }
}
