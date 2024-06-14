package com.hartwig.hmftools.bamtools.metrics;

public enum FlagStatType
{
    TOTAL,
    PRIMARY, // not a supplementary or secondary read
    SECONDARY,
    SUPPLEMENTARY,
    DUPLICATE, // includes supplementaries
    PRIMARY_DUPLICATE,
    MAPPED,
    PRIMARY_MAPPED,
    PAIRED,
    READ1,
    READ2,
    PROPERLY_PAIRED,
    PAIR_MAPPED,
    SINGLETON,
    INTER_CHR_PAIR_MAPPED,
    INTER_CHR_PAIR_MAP_QUAL_GE5;
}
