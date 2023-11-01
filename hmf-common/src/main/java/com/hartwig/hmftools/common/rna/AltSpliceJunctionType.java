package com.hartwig.hmftools.common.rna;

public enum AltSpliceJunctionType
{
    SKIPPED_EXONS,
    NOVEL_5_PRIME,
    NOVEL_3_PRIME,
    NOVEL_EXON,
    NOVEL_INTRON,
    MIXED_TRANS,
    INTRONIC,
    EXON_INTRON,
    CIRCULAR,
    CANONICAL, // used to put known splice sites into the same form as alt-SJs for Cuppa
    UNKNOWN
}
