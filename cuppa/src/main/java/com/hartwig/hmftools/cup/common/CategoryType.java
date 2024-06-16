package com.hartwig.hmftools.cup.common;

import java.util.Arrays;
import java.util.List;

public enum CategoryType
{
    SNV,
    SV,
    SAMPLE_TRAIT,
    FEATURE,
    GENE_EXP,
    ALT_SJ,

    @Deprecated
    COMBINED;

    @Deprecated
    public static boolean isDna(final CategoryType type)
    {
        return type == SNV || type == SV || type == SAMPLE_TRAIT || type == FEATURE;
    }

    @Deprecated
    public static boolean isRna(final CategoryType type)
    {
        return type == GENE_EXP || type == ALT_SJ;
    }

    public static List<CategoryType> getDnaCategories()
    {
        return Arrays.asList(SNV, SV, SAMPLE_TRAIT, FEATURE);
    }

    public static List<CategoryType> getRnaCategories()
    {
        return Arrays.asList(GENE_EXP, ALT_SJ);
    }

    public static List<CategoryType> getAllCategories()
    {
        return Arrays.asList(SNV, SV, SAMPLE_TRAIT, FEATURE, GENE_EXP, ALT_SJ);
    }
}
