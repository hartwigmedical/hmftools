package com.hartwig.hmftools.cup.common;

public enum ClassifierType
{
    SNV_96_PAIRWISE_SIMILARITY,
    GENOMIC_POSITION_SIMILARITY,
    EXPRESSION_PAIRWISE,
    EXPRESSION_COHORT,
    FEATURE,
    ALT_SJ_COHORT,
    ALT_SJ_PAIRWISE,
    COMBINED;

    public static boolean isDna(final ClassifierType type)
    {
        return type == SNV_96_PAIRWISE_SIMILARITY || type == GENOMIC_POSITION_SIMILARITY || type == FEATURE;
    }

    public static boolean isRna(final ClassifierType type)
    {
        return type == EXPRESSION_COHORT || type == EXPRESSION_PAIRWISE || type == ALT_SJ_COHORT || type == ALT_SJ_PAIRWISE;
    }
}
