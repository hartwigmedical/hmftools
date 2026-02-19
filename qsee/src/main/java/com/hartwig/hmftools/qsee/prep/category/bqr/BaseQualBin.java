package com.hartwig.hmftools.qsee.prep.category.bqr;

import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.sequencing.SequencingType;

public enum BaseQualBin
{
    LOW,
    MEDIUM,
    HIGH;

    public boolean withinBin(byte baseQual, SequencingType sequencingType)
    {
        return switch(this)
        {
            case LOW -> BaseQualAdjustment.isLowBaseQual(baseQual);
            case MEDIUM -> BaseQualAdjustment.isMediumBaseQual(baseQual, sequencingType);
            case HIGH -> BaseQualAdjustment.isHighBaseQual(baseQual, sequencingType);
        };
    }
}


