package com.hartwig.hmftools.cup.common;

import com.hartwig.hmftools.common.purple.Gender;

public class CupConstants
{
    public static final double SNV_CSS_THRESHOLD = 0.01;

    public static final double CSS_SIMILARITY_CUTOFF = 0.01;
    public static final int CSS_SIMILARITY_MAX_MATCHES = 20;

    public static final double SNV_POS_FREQ_CSS_THRESHOLD = 0.01;
    public static final double RNA_GENE_EXP_CSS_THRESHOLD = 0.01;

    public static final double SNV_CSS_DIFF_EXPONENT = 8;
    public static final double SNV_POS_FREQ_DIFF_EXPONENT = 10;
    public static final double RNA_GENE_EXP_DIFF_EXPONENT = 50;
    public static final double RNA_ALT_SJ_DIFF_EXPONENT = 3.5;

    public static final int POS_FREQ_BUCKET_SIZE = 500000;
    public static final int POS_FREQ_MAX_SAMPLE_COUNT = 20000;

    public static final double MIN_CLASSIFIER_SCORE = 0.01;

    public static final double FEATURE_DAMPEN_FACTOR = 0.8;
    public static final double COMBINED_DAMPEN_FACTOR = 0.4;
    public static final double DNA_DAMPEN_FACTOR = 0.65;
    public static final double RNA_DAMPEN_FACTOR = 0.7;

    public static final String CANCER_TYPE_UNKNOWN = "Unknown";
    public static final String CANCER_TYPE_OTHER = "Other";
    public static final String CANCER_TYPE_PAN = "ALL";

    public static final double DRIVER_ZERO_PREVALENCE_ALLOCATION = 0.10;
    public static final double NON_DRIVER_ZERO_PREVALENCE_ALLOCATION = 0.01;

    // cancer types with gender-exclusions
    public static final String CANCER_TYPE_PROSTATE = "Prostate";
    public static final String CANCER_TYPE_OVARY = "Ovary";
    public static final String CANCER_TYPE_UTERUS = "Uterus";
    public static final String CANCER_TYPE_TESTIS = "Testis";

    public static boolean isCandidateCancerType(final Gender gender, final String cancerType)
    {
        if(gender == null)
            return true;

        if(cancerType.contains(CANCER_TYPE_UTERUS) || cancerType.contains(CANCER_TYPE_OVARY))
        {
            if(gender != Gender.FEMALE)
                return false;
        }

        if(cancerType.contains(CANCER_TYPE_PROSTATE) || cancerType.contains(CANCER_TYPE_TESTIS))
        {
            if(gender == Gender.FEMALE)
                return false;
        }

        return true;
    }
}
