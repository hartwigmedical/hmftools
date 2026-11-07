package com.hartwig.hmftools.wisp.variant;

public enum FilterReason
{
    NO_PASS,
    GERMLINE_AF,
    NEARBY_INDEL,
    NON_SNV,
    MAPPABILITY,
    REPEAT_COUNT,
    LOW_CONFIDENCE,
    SUBCLONAL,
    GC_RATIO,
    LOW_QUAL_PER_AD,
    AVG_EDGE_DIST,
    OUTLIER,
    DUAL_ERROR_RATE;

    public static boolean isMosiacType(final FilterReason filterReason)
    {
        switch(filterReason)
        {
            case NO_PASS:
            case MAPPABILITY:
            case REPEAT_COUNT:
            case LOW_CONFIDENCE:
            case GC_RATIO:
            case OUTLIER:
            case DUAL_ERROR_RATE:
                return true;

            default:
                return false;
        }
    }
}
