package com.hartwig.hmftools.cup.common;

public enum ClassifierType
{
    SNV_96_PAIRWISE_SIMILARITY,
    GENOMIC_POSITION_SIMILARITY,
    EXPRESSION_PAIRWISE,
    EXPRESSION_COHORT,
    FEATURE,
    COMBINED;

    public static boolean isDna(final ClassifierType type)
    {
        return type == SNV_96_PAIRWISE_SIMILARITY || type == GENOMIC_POSITION_SIMILARITY || type == FEATURE;
    }

    public static boolean isRna(final ClassifierType type)
    {
        return type == EXPRESSION_COHORT || type == EXPRESSION_PAIRWISE;
    }
}
