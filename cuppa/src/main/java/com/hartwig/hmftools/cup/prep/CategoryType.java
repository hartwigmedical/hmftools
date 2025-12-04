package com.hartwig.hmftools.cup.prep;

import java.util.Arrays;
import java.util.List;

public enum CategoryType
{
    SNV(SourceDataType.DNA),
    SV(SourceDataType.DNA),
    SAMPLE_TRAIT(SourceDataType.DNA),
    DRIVER(SourceDataType.DNA),

    GENE_EXP(SourceDataType.RNA),
    ALT_SJ(SourceDataType.RNA);

    final SourceDataType mSourceDataType;

    CategoryType(SourceDataType sourceDataType)
    {
        mSourceDataType = sourceDataType;
    }

    public static List<CategoryType> getDnaCategories()
    {
        return Arrays.stream(CategoryType.values()).filter(x -> x.mSourceDataType == SourceDataType.DNA).toList();
    }

    public static List<CategoryType> getRnaCategories()
    {
        return Arrays.stream(CategoryType.values()).filter(x -> x.mSourceDataType == SourceDataType.RNA).toList();
    }

    public static List<CategoryType> getAllCategories()
    {
        return List.of(CategoryType.values());
    }

    public boolean isDna() { return this.mSourceDataType == SourceDataType.DNA; }
    public boolean isRna() { return this.mSourceDataType == SourceDataType.RNA; }

    private enum SourceDataType
    {
        DNA,
        RNA
    }
}
