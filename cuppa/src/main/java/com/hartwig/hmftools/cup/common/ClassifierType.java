package com.hartwig.hmftools.cup.common;

public enum ClassifierType
{
    SNV_96_PAIRWISE_SIMILARITY,
    GENOMIC_POSITION_SIMILARITY,
    GENE_EXPRESSION_PAIRWISE,
    GENE_EXPRESSION_CANCER,
    FEATURE_PREVALENCE,
    COMBINED;

    public static boolean isDna(final ClassifierType type)
    {
        return type == SNV_96_PAIRWISE_SIMILARITY || type == GENOMIC_POSITION_SIMILARITY || type == FEATURE_PREVALENCE;
    }

    public static boolean isRna(final ClassifierType type)
    {
        return type == GENE_EXPRESSION_CANCER || type == GENE_EXPRESSION_PAIRWISE;
    }
}
