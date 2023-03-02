package com.hartwig.hmftools.datamodel.purple;

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
}
