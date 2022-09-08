package com.hartwig.hmftools.compar;

import java.util.List;

import com.google.common.collect.Lists;

public enum Category
{
    PURITY,
    DRIVER,
    COPY_NUMBER,
    GENE_COPY_NUMBER,
    SOMATIC_VARIANT,
    GERMLINE_VARIANT,
    GERMLINE_DELETION,
    FUSION,
    DISRUPTION,
    GERMLINE_SV,
    CUPPA,
    LILAC,
    CHORD;

    public static final String ALL_CATEGORIES = "ALL";
    public static final String LINX_CATEGORIES = "LINX";
    public static final String PURPLE_CATEGORIES = "PURPLE";

    public static List<Category> purpleCategories()
    {
        return Lists.newArrayList(PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT ,GERMLINE_DELETION);
    }

    public static List<Category> linxCategories() { return Lists.newArrayList(DRIVER, FUSION, DISRUPTION, GERMLINE_SV); }
}
