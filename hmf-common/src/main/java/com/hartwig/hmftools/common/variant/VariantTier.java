package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public enum VariantTier
{
    UNKNOWN,
    LOW_CONFIDENCE,
    HIGH_CONFIDENCE,
    PANEL,
    HOTSPOT;

    public static final String TIER = "TIER";

    @NotNull
    public static VariantTier fromContext(final VariantContext context)
    {
        return fromString(context.getAttributeAsString(TIER, UNKNOWN.toString()));
    }

    @NotNull
    public static VariantTier fromString(final String string)
    {
        switch(string)
        {
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
