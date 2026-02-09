package com.hartwig.hmftools.isofox;

import java.util.List;

public enum IsofoxFunction
{
    TRANSCRIPT_COUNTS,
    ALT_SPLICE_JUNCTIONS,
    RETAINED_INTRONS,
    FUSIONS,
    STATISTICS,
    READ_COUNTS,
    NEO_EPITOPES,
    UNMAPPED_READS;

    public static final List<IsofoxFunction> DEFAULT_FUNCTIONS = List.of(TRANSCRIPT_COUNTS, ALT_SPLICE_JUNCTIONS, FUSIONS);

}
