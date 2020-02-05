package com.hartwig.hmftools.svtools.rna_expression;

public enum RegionMatchType
{
    NONE,
    EXON_BOUNDARY,  // read matches one exon boundary
    WITHIN_EXON,    // read fully contained within the exon
    EXON_MATCH,     // read fully contained within the exon
    EXON_INTRON,      // reads spanning to unmapped regions where adjacent regions exist
    INTRONIC;
}
