package com.hartwig.hmftools.cup.common;

public class CupConstants
{
    public static final double SNV_CSS_THRESHOLD = 0.01;

    public static final double SNV_CSS_SIMILARITY_CUTOFF = 0.8;
    public static final int SNV_CSS_SIMILARITY_MAX_MATCHES = 5;

    public static final double SNV_POS_FREQ_CSS_THRESHOLD = 0.01;
    public static final double RNA_GENE_EXP_CSS_THRESHOLD = 0.01;

    public static final double SNV_CSS_DIFF_EXPONENT = 8;
    public static final double SNV_POS_FREQ_DIFF_EXPONENT = 10;
    public static final double RNA_GENE_EXP_DIFF_EXPONENT = 10;

    public static final int POS_FREQ_BUCKET_SIZE = 500000;
    public static final int POS_FREQ_MAX_SAMPLE_COUNT = 20000;


    public static final double MIN_CLASSIFIER_SCORE = 0.02;

    public static final String CANCER_TYPE_UNKNOWN = "Unknown";

    // cancer types with gender-exclusions
    public static final String CANCER_TYPE_PROSTATE = "Prostate";
    public static final String CANCER_TYPE_OVARY = "Ovary";
    public static final String CANCER_TYPE_UTERUS = "Uterus";
    public static final String CANCER_TYPE_PAN = "ALL";

    public static final double DRIVER_ZERO_PREVALENCE_ALLOCATION = 0.10;
    public static final double NON_DRIVER_ZERO_PREVALENCE_ALLOCATION = 0.02;

}
