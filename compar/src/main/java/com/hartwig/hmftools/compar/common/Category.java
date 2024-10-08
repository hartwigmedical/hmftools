package com.hartwig.hmftools.compar.common;

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
    CHORD,
    PEACH,
    VIRUS,
    TUMOR_FLAGSTAT,
    GERMLINE_FLAGSTAT,
    TUMOR_BAM_METRICS,
    GERMLINE_BAM_METRICS,
    SNP_GENOTYPE;

    public static final String ALL_CATEGORIES = "ALL";
    public static final String LINX_CATEGORIES = "LINX";
    public static final String PURPLE_CATEGORIES = "PURPLE";
    public static final String PANEL_CATEGORIES = "PANEL";

    public static List<Category> purpleCategories()
    {
        return Lists.newArrayList(PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT ,GERMLINE_DELETION);
    }

    public static List<Category> linxCategories() { return Lists.newArrayList(DRIVER, FUSION, DISRUPTION, GERMLINE_SV); }

    public static List<Category> panelCategories()
    {
        return Lists.newArrayList(PURITY, DRIVER, SOMATIC_VARIANT, FUSION, DISRUPTION, TUMOR_FLAGSTAT, TUMOR_BAM_METRICS, SNP_GENOTYPE);
    }
}
