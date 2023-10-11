package com.hartwig.hmftools.wisp.purity;

public class PurityConstants
{
    public static final double MAX_SUBCLONAL_LIKELIHOOD = 0.1;
    public static final double MAX_REPEAT_COUNT = 3;

    public static final int MIN_QUAL_PER_AD = 18;
    public static final int LOW_QUAL_NOISE_CUTOFF = 25;

    public static final double MAX_COPY_NUMBER = 6;
    public static final double CLONAL_COPY_NUMBER_MARGIN = 0.2;

    public static final double DEFAULT_GC_RATIO_MIN = 0.4;

    public static final double DEFAULT_NOISE_READS_PER_MILLION = 30;
    public static final double DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND = 1;

    public static final double SOMATIC_PEAK_MAX_PROBABILITY = 0.05;
    public static final int SOMATIC_PEAK_MIN_VARIANTS = 10;
    public static final int SOMATIC_PEAK_MIN_DEPTH = 100;
    public static final int SOMATIC_PEAK_MIN_AD = 5;
    public static final int SOMATIC_PEAK_MIN_PEAK_VARIANTS = 5;

    public static final int VAF_PEAK_MODEL_MIN_FRAG_VARIANTS = 10;
    public static final int VAF_PEAK_MODEL_MIN_AVG_DEPTH = 20;
    public static final int LOW_COUNT_MODEL_MIN_FRAG_VARIANTS = 5;
    public static final int LOW_COUNT_MODEL_MIN_AVG_DEPTH = 50;
    public static final double DROPOUT_RATE_INCREMENT = 0.1;

    public static final String PURPLE_CTDNA_SOMATIC_VCF_ID = ".purple.somatic.ctdna.";
}
