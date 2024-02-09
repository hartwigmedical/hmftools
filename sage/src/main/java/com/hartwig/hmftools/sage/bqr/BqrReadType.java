package com.hartwig.hmftools.sage.bqr;

import com.hartwig.hmftools.common.samtools.UmiReadType;

public enum BqrReadType
{
    NONE,
    SINGLE,
    DUAL,
    ULT_BALANCED,
    ULT_STANDARD;

    public static BqrReadType fromUmiType(final UmiReadType umiReadType)
    {
        switch(umiReadType)
        {
            case NONE: return BqrReadType.NONE;
            case SINGLE: return BqrReadType.SINGLE;
            case DUAL: return BqrReadType.DUAL;
            default: return BqrReadType.NONE;
        }
    }
}
