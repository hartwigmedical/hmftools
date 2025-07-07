package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE;

import javax.annotation.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public enum StructuralVariantType
{
    BND,
    DEL,
    DUP,
    INF,
    INS,
    INV,
    SGL;

    public static @Nullable StructuralVariantType fromContext(final VariantContext context)
    {
        return context.hasAttribute(SV_TYPE) ? StructuralVariantType.valueOf(context.getAttributeAsString(SV_TYPE, "")) : null;
    }

    public static int typeAsInt(StructuralVariantType type) { return type.ordinal(); }
}
