package com.hartwig.hmftools.cup.common;

public enum ClassifierType
{
    SNV_COUNT_CSS,
    SNV_POS_FREQ_CSS,
    FEATURE_PREVALENCE,
    PERCENTILES;

    public static String displayString(final ClassifierType type)
    {
        switch(type)
        {
            case SNV_COUNT_CSS: return "SNV_COUNT_CSS";
            case SNV_POS_FREQ_CSS: return "GEN_POS_CSS";
            case FEATURE_PREVALENCE: return "FEATURE_PREV";
            case PERCENTILES: return "PERCENTILES";
            default: return "UNKNOWN";
        }
    }

}
