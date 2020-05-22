package com.hartwig.hmftools.linx.fusion.rna;

public enum RnaJunctionType
{
    NOT_SET,
    KNOWN, // matches known exon boundary (ie from Ensembl)
    CANONICAL, // matches known motif
    UNKNOWN;
}
