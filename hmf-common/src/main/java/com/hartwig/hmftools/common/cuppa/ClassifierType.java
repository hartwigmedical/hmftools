package com.hartwig.hmftools.common.cuppa;

import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_COMBINED;
import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_DNA_COMBINED;
import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_RNA_COMBINED;

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
        if(classiferType.equals(DATA_TYPE_DNA_COMBINED) || classiferType.equals(DATA_TYPE_RNA_COMBINED) || classiferType.equals(DATA_TYPE_COMBINED))
            return COMBINED;

        return ClassifierType.valueOf(classiferType);
    }

}
