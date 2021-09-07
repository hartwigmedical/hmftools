package com.hartwig.hmftools.common.sv;

import org.jetbrains.annotations.NotNull;

public enum StructuralVariantType {
    BND,
    DEL,
    DUP,
    INF,
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

    public static int typeAsInt(@NotNull StructuralVariantType type)
    {
        // ordered alphabetically
        switch(type)
        {
            case BND: return 0;
            case DEL: return 1;
            case DUP: return 2;
            case INF: return 3;
            case INS: return 4;
            case INV: return 5;
            case SGL: return 6;
        }

        return 0;
    }
}
