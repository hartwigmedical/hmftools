package com.hartwig.hmftools.sage.bqr;

import com.hartwig.hmftools.common.samtools.UmiReadType;
import com.hartwig.hmftools.common.sequencing.UltimaConsensusType;

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

    public static BqrReadType fromUltimaType(final UltimaConsensusType consensusType)
    {
        switch(consensusType)
        {
            case STANDARD: return BqrReadType.NONE;
            case BALANCED: return BqrReadType.DUAL;
            default: return BqrReadType.NONE;
        }
    }
}
