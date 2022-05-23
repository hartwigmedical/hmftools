package com.hartwig.hmftools.cup.common;

public enum ClassifierType
{
    SNV_96_PAIRWISE,
    GENOMIC_POSITION_COHORT,
    GENOMIC_POSITION_PAIRWISE,
    EXPRESSION_PAIRWISE,
    EXPRESSION_COHORT,
    FEATURE,
    GENDER,
    ALT_SJ_COHORT,
    ALT_SJ_PAIRWISE;

    public static boolean isDna(final ClassifierType type)
    {
        return type == SNV_96_PAIRWISE || type == GENOMIC_POSITION_COHORT || type == GENOMIC_POSITION_PAIRWISE || type == FEATURE || type == GENDER;
    }

    public static boolean isRna(final ClassifierType type)
    {
        return type == EXPRESSION_COHORT || type == EXPRESSION_PAIRWISE || type == ALT_SJ_COHORT || type == ALT_SJ_PAIRWISE;
    }

    public static boolean applyMinScore(final ClassifierType type)
    {
        return type != GENDER;
    }

}
