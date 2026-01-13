package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public enum HotspotType
{
    HOTSPOT,
    NEAR_HOTSPOT,
    NON_HOTSPOT;

    public static final String HOTSPOT_FLAG = HOTSPOT.toString();
    public static final String NEAR_HOTSPOT_FLAG = NEAR_HOTSPOT.toString();

    public static final String HOTSPOT_DESCRIPTION = "Site is at a known hotspot location";

    @NotNull
    public static HotspotType fromVariant(@NotNull final VariantContext context)
    {
        if(context.getAttributeAsBoolean(HOTSPOT_FLAG, false))
        {
            return HotspotType.HOTSPOT;
        }

        if(context.getAttributeAsBoolean(NEAR_HOTSPOT_FLAG, false))
        {
            return HotspotType.NEAR_HOTSPOT;
        }

        return HotspotType.NON_HOTSPOT;
    }

}
