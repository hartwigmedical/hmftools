package com.hartwig.hmftools.common.purple.segment;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

public enum StructuralVariantSupport {
    BND,
    INV,
    DEL,
    INS,
    DUP,
    MULTIPLE,
    NONE;

    public static StructuralVariantSupport fromVariant(StructuralVariantType type) {
        return StructuralVariantSupport.valueOf(type.toString());
    }
}
