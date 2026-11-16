package com.hartwig.hmftools.isofox;

import java.util.List;

public enum IsofoxFunction
{
    TRANSCRIPT_COUNTS,
    ALT_SPLICE_JUNCTIONS,
    RETAINED_INTRONS,
    FUSIONS,
    STATISTICS,
    NEO_EPITOPES;

    public static final List<IsofoxFunction> DEFAULT_FUNCTIONS = List.of(TRANSCRIPT_COUNTS, ALT_SPLICE_JUNCTIONS, FUSIONS);

}
