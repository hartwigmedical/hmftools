package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.pow;

public class PurityConstants
{
    public static final double MAX_SUBCLONAL_LIKELIHOOD = 0.1;
    public static final double MAX_REPEAT_COUNT = 3;

    public static final int MIN_QUAL_PER_AD = 18;
    public static final int LOW_QUAL_NOISE_CUTOFF = 25;

    public static final double MAX_COPY_NUMBER = 6;
    public static final double CLONAL_COPY_NUMBER_MARGIN = 0.2;

    public static final double VARIANT_OUTLIER_VAF_MULTIPLE = 8;
    public static final int VARIANT_OUTLIER_MIN_AD = 5;
    public static final double VARIANT_OUTLIER_MIN_AD_PERC = 0.1;

    public static final double DEFAULT_GC_RATIO_MIN = 0.4;

    public static final double DEFAULT_NOISE_READS_PER_MILLION = 30;
    public static final double DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND = 1;

    public static final double DROPOUT_RATE_VAF_INCREMENT = 0.05;
    public static final int DROPOUT_RATE_MIN_DEPTH = 100;
    public static final double DROPOUT_RATE_PROBABILITY = pow(10, -5);

    public static final String PURPLE_CTDNA_SOMATIC_VCF_ID = ".purple.somatic.ctdna.";
}
