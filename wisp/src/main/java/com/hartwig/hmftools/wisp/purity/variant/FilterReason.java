package com.hartwig.hmftools.wisp.purity.variant;

public enum FilterReason
{
    NO_FILTER,
    NO_PASS,
    NON_SNV,
    MAPPABILITY,
    REPEAT_COUNT,
    LOW_CONFIDENCE,
    SUBCLONAL,
    GC_RATIO,
    LOW_QUAL_PER_AD,
    CHIP;
}
