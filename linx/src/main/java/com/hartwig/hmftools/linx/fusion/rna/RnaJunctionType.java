package com.hartwig.hmftools.linx.fusion.rna;

public enum RnaJunctionType
{
    NOT_SET,
    KNOWN, // matches known exon boundary (ie from Ensembl)
    CANONICAL, // matches known motif
    UNKNOWN;

    public static boolean isUnspliced(RnaJunctionType type) { return type == UNKNOWN || type == NOT_SET; }
}
