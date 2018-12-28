package com.hartwig.hmftools.common.variant.structural;

import org.jetbrains.annotations.NotNull;

public enum StructuralVariantType {
    BND,
    DEL,
    DUP,
    INS,
    INV,
    SGL;

    @NotNull
    public static StructuralVariantType fromAttribute(@NotNull String svType) {
        if (svType.startsWith("DUP")) {
            return DUP;
        }
        return StructuralVariantType.valueOf(svType);
    }

    public static int SV_TYPE_COUNT = 6;

    public static int typeAsInt(StructuralVariantType type)
    {
        // ordered alphabetically
        switch(type)
        {
            case BND: return 0;
            case DEL: return 1;
            case DUP: return 2;
            case INS: return 3;
            case INV: return 4;
            case SGL: return 5;
        }

        return 0;
    }
}
