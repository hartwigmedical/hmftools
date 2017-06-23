package com.hartwig.hmftools.common.purple.segment;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

public enum PurpleSegmentSource {
    BND,
    INV,
    DEL,
    INS,
    DUP,
    FREEC;

    public static PurpleSegmentSource fromVariant(StructuralVariantType type) {
        return PurpleSegmentSource.valueOf(type.toString());
    }
}
