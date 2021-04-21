package com.hartwig.hmftools.linx.types;

public enum ResolvedType
{
    NONE,

    // artefacts
    DUP_BE,
    LOW_VAF,
    PAIR_INF,
    SGL_MAPPED_INF,

    // simple types including synthetic instances
    DEL,
    DUP,
    INF,
    SGL,
    INV,
    INS,
    UNBAL_TRANS,
    SIMPLE_GRP, // group of DUPs, INS and DELs - always dissolved in final classification

    // 2-SV types:
    RECIP_TRANS,
    RECIP_TRANS_DEL_DUP,
    RECIP_TRANS_DUPS,
    RECIP_INV,
    RECIP_INV_DUPS,
    RECIP_INV_DEL_DUP,
    DUP_TI,
    DEL_TI,
    UNBAL_TRANS_TI,
    PAIR_OTHER,
    RESOLVED_FOLDBACK,
    FB_INV_PAIR,
    SGL_PAIR_DEL,
    SGL_PAIR_DUP,
    SGL_PAIR_INS,

    // complex types
    LINE,
    COMPLEX,
    DOUBLE_MINUTE;


    public boolean isSimple()
    {
        return (this == ResolvedType.DEL || this == ResolvedType.DUP || this == ResolvedType.INS
                || this == SGL_PAIR_DEL || this == SGL_PAIR_DUP || this == SGL_PAIR_INS);
    }

}
