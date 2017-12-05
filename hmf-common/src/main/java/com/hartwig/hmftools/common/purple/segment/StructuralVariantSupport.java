package com.hartwig.hmftools.common.purple.segment;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

public enum StructuralVariantSupport {
    BND,
    INV,
    DEL,
    INS,
    DUP,
    MULTIPLE,
    NONE;

    @NotNull
    public static StructuralVariantSupport fromVariant(@NotNull StructuralVariantType type) {
        return StructuralVariantSupport.valueOf(type.toString());
    }
}
