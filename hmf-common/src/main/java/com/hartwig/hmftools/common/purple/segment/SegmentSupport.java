package com.hartwig.hmftools.common.purple.segment;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

public enum SegmentSupport {
    BND(true),
    INV(true),
    DEL(true),
    SGL(true),
    INS(true),
    INF(true),
    DUP(true),
    MULTIPLE(true),
    CENTROMERE(false),
    TELOMERE(false),
    NONE(false),
    UNKNOWN(false);

    private final boolean isSV;

    SegmentSupport(boolean isSV) {
        this.isSV = isSV;
    }

    public boolean isSV() {
        return isSV;
    }

    public static SegmentSupport fromVariant(StructuralVariantType type) {
        switch (type) {
            case INF: return INF;
            case BND: return BND;
            case DEL: return DEL;
            case INS: return INS;
            case DUP: return DUP;
            case INV: return INV;
            case SGL: return SGL;
        }

        throw new IllegalArgumentException("Unknown variant type: " + type);
    }

}
