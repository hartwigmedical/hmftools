package com.hartwig.hmftools.common.qual;

import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.common.sequencing.UltimaConsensusType;

public enum BqrReadType
{
    NONE(false),
    SINGLE(false),
    DUAL(true),
    ULT_BALANCED(true),
    ULT_STANDARD(false);

    private final boolean mIsHighQuality;

    BqrReadType(boolean isHighQuality) { mIsHighQuality = isHighQuality; }

    public boolean isHighQuality() { return mIsHighQuality; }

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
