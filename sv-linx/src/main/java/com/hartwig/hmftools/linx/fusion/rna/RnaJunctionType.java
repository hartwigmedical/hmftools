package com.hartwig.hmftools.linx.fusion.rna;

public enum RnaJunctionType
{
    KNOWN, // matches known exon boundary (ie from Ensembl)
    CANONICAL, // matches known motif
    UNKNOWN;
}
