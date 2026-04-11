package com.hartwig.hmftools.linx.cohort;

import com.hartwig.hmftools.linx.types.ResolvedType;

public enum SvType
{
    INS,
    DEL_100,
    DEL_100_1K,
    DEL_1_10K,
    DEL_10_100K,
    DEL_LONG,
    DUP_100,
    DUP_100_1K,
    DUP_1_10K,
    DUP_10_100K,
    DUP_LONG,
    BND,
    INV,
    SGL,

    LINE,

    // clustered types
    DEL_SYNTHETIC,
    DUP_SYNTHETIC,
    RECIP_TRANS,
    RECIP_INV,
    PAIR_OTHER,
    COMPLEX;

    public static boolean ignoreResolvedType(final ResolvedType resolvedType)
    {
        switch(resolvedType)
        {
            case DUP_BE:
            case INF:
            case LOW_VAF:
                return true;

            default:
                return false;
        }
    }

    public static SvType classifySv(final SvData sv)
    {
        if(sv.ClusterCount == 1)
        {
            switch(sv.SvType)
            {
                case DEL:
                    if(sv.length() < 100)
                        return DEL_100;
                    else if(sv.length() < 1_000)
                        return DEL_100_1K;
                    else if(sv.length() < 10_000)
                        return DEL_1_10K;
                    else if(sv.length() < 100_000)
                        return DEL_10_100K;
                    else
                        return DEL_LONG;

                case DUP:
                    if(sv.length() < 100)
                        return DUP_100;
                    else if(sv.length() < 1_000)
                        return DUP_100_1K;
                    else if(sv.length() < 10_000)
                        return DUP_1_10K;
                    else if(sv.length() < 100_000)
                        return DUP_10_100K;
                    else
                        return DUP_LONG;

                case INS: return SvType.INS;
                case INV: return SvType.INV;
                case BND: return SvType.BND;

                case SGL:
                default:
                    return SvType.SGL;
            }
        }

        switch(sv.Resolved)
        {
            case RECIP_TRANS:
            case RECIP_TRANS_DUPS:
            case RECIP_TRANS_DEL_DUP:
                return RECIP_TRANS;

            case RECIP_INV:
            case RECIP_INV_DUPS:
            case RECIP_INV_DEL_DUP:
                return RECIP_INV;

            case DUP_TI:
            case SGL_PAIR_DUP:
            case SGL_PAIR_INS:
                return DUP_SYNTHETIC;

            case DEL_TI:
            case SGL_PAIR_DEL:
                return DEL_SYNTHETIC;

            case LINE:
                return LINE;

            default:
                if(sv.ClusterCount == 2)
                    return PAIR_OTHER;
                else
                    return COMPLEX;
        }
    }
}
