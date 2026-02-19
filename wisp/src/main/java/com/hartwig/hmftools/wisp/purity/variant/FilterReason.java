package com.hartwig.hmftools.wisp.purity.variant;

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
    OUTLIER;
}
