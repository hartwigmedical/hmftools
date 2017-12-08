package com.hartwig.hmftools.common.variant.structural;

import org.jetbrains.annotations.NotNull;

public enum StructuralVariantType {
    BND,
    INV,
    DEL,
    INS,
    DUP;

    @NotNull
    public static StructuralVariantType fromAttribute(@NotNull String svType) {
        if (svType.startsWith("DUP")) {
            return DUP;
        }
        return StructuralVariantType.valueOf(svType);
    }
}
