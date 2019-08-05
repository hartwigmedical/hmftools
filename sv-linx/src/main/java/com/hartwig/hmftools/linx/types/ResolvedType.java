package com.hartwig.hmftools.linx.types;

public enum ResolvedType
{
    NONE,
    DUP_BE,
    LOW_VAF,
    LINE,
    DEL,
    DUP,
    INF,
    SGL,
    INV,
    INS,
    RECIP_TRANS,
    RECIP_INV,
    UNBAL_TRANS,
    RECIP_DUPS,
    RECIP_DUP_DEL,
    COMPLEX,
    PAIR_OTHER,
    SIMPLE_GRP,
    FB_INV_PAIR,
    SGL_PAIR_DEL,
    SGL_PAIR_DUP,
    SGL_PAIR_INS;

    public boolean isSimple()
    {
        return (this == ResolvedType.DEL || this == ResolvedType.DUP || this == ResolvedType.INS
                || this == SGL_PAIR_DEL || this == SGL_PAIR_DUP || this == SGL_PAIR_INS);
    }

}
