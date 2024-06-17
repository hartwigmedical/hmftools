package com.hartwig.hmftools.cup.prep;

import java.util.List;

public enum CategoryType
{
    SNV,
    SV,
    SAMPLE_TRAIT,
    DRIVER,
    GENE_EXP,
    ALT_SJ;

    public static List<CategoryType> getDnaCategories()
    {
        return List.of(SNV, SV, SAMPLE_TRAIT, DRIVER);
    }

    public static List<CategoryType> getRnaCategories()
    {
        return List.of(GENE_EXP, ALT_SJ);
    }

    public static List<CategoryType> getAllCategories()
    {
        return List.of(SNV, SV, SAMPLE_TRAIT, DRIVER, GENE_EXP, ALT_SJ);
    }
}
