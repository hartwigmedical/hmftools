package com.hartwig.hmftools.common.cuppa;

public enum ClassifierType
{
    SNV_96_PAIRWISE,
    GENOMIC_POSITION_COHORT,
    EXPRESSION_PAIRWISE,
    EXPRESSION_COHORT,
    FEATURE,
    GENDER,
    ALT_SJ_COHORT,
    ALT_SJ_PAIRWISE,
    COMBINED;

    public static boolean isDna(final ClassifierType type)
    {
        return type == SNV_96_PAIRWISE || type == GENOMIC_POSITION_COHORT || type == FEATURE || type == GENDER;
    }

    public static boolean isRna(final ClassifierType type)
    {
        return type == EXPRESSION_COHORT || type == EXPRESSION_PAIRWISE || type == ALT_SJ_COHORT || type == ALT_SJ_PAIRWISE;
    }

    public static boolean applyMinScore(final ClassifierType type)
    {
        return type != GENDER;
    }

    public static ClassifierType fromString(final String classiferType)
    {
        if(classiferType.equals("SNV_96_PAIRWISE_SIMILARITY")) // backwards compatibility (earlier than 1.7)
            return SNV_96_PAIRWISE;

        if(classiferType.equals("GENOMIC_POSITION_SIMILARITY")) // as above
            return GENOMIC_POSITION_COHORT;

        return ClassifierType.valueOf(classiferType);
    }

}
