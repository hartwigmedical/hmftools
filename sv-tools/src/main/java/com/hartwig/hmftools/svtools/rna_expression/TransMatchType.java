package com.hartwig.hmftools.svtools.rna_expression;

// type of match against a transcript
public enum TransMatchType
{
    UNKNOWN,
    SPLICE_JUNCTION,
    EXONIC,
    UNSPLICED,
    OTHER_TRANS,
    ALT;
}
