package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public enum Hotspot
{
    HOTSPOT,
    NEAR_HOTSPOT,
    NON_HOTSPOT;

    public static final String HOTSPOT_FLAG = HOTSPOT.toString();
    public static final String NEAR_HOTSPOT_FLAG = NEAR_HOTSPOT.toString();

    @NotNull
    public static Hotspot fromVariant(@NotNull final VariantContext context)
    {
        if(context.getAttributeAsBoolean(HOTSPOT_FLAG, false))
        {
            return Hotspot.HOTSPOT;
        }

        if(context.getAttributeAsBoolean(NEAR_HOTSPOT_FLAG, false))
        {
            return Hotspot.NEAR_HOTSPOT;
        }

        return Hotspot.NON_HOTSPOT;
    }

}
