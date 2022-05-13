package com.hartwig.hmftools.cup.common;

public enum CategoryType
{
    SNV,
    SV,
    SAMPLE_TRAIT,
    FEATURE,
    GENE_EXP,
    ALT_SJ,
    CLASSIFIER,
    COMBINED;

    public static final String ALL_CATEGORIES = "ALL";
    public static final String DNA_CATEGORIES = "DNA";
    public static final String RNA_CATEGORIES = "RNA";

    public static boolean isDna(final CategoryType type)
    {
        return type == SNV || type == SV || type == SAMPLE_TRAIT || type == FEATURE;
    }

    public static boolean isRna(final CategoryType type)
    {
        return type == GENE_EXP || type == ALT_SJ;
    }

    public static boolean isSummary(final CategoryType type) { return type == CLASSIFIER || type == COMBINED; }

}
