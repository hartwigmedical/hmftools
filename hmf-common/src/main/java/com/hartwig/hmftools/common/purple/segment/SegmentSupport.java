package com.hartwig.hmftools.common.purple.segment;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

public enum SegmentSupport {
    BND,
    INV,
    DEL,
    INS,
    DUP,
    MULTIPLE,
    CENTROMERE,
    TELOMERE,
    NONE;

    public static SegmentSupport fromVariant(StructuralVariantType type) {
        return SegmentSupport.valueOf(type.toString());
    }

    public static SegmentSupport fromSVSupport(StructuralVariantSupport support) {
        switch (support) {
            case BND:
                return BND;
            case INV:
                return INV;
            case DEL:
                return DEL;
            case INS:
                return INS;
            case DUP:
                return DUP;
            case MULTIPLE:
                return MULTIPLE;
            default:
                return NONE;
        }
    }

}
