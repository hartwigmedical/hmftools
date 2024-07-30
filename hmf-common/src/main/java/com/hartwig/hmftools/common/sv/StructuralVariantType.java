package com.hartwig.hmftools.common.sv;

import org.jetbrains.annotations.NotNull;

public enum StructuralVariantType
{
    BND,
    DEL,
    DUP,
    INF,
    INS,
    INV,
    SGL;

    @NotNull
    public static StructuralVariantType fromAttribute(@NotNull String svType)
    {
        if(svType.startsWith("DUP"))
        {
            return DUP;
        }
        return StructuralVariantType.valueOf(svType);
    }

    public static int typeAsInt(StructuralVariantType type) { return type.ordinal(); }
}
