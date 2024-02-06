package com.hartwig.hmftools.wisp.purity;

public class PurityConstants
{
    public static final double MAX_SUBCLONAL_LIKELIHOOD = 0.5;
    public static final double SUBCLONAL_VCN_THRESHOLD = 0.7;
    public static final double MAX_REPEAT_COUNT = 3;

    public static final int MIN_QUAL_PER_AD = 18;
    public static final int LOW_QUAL_NOISE_CUTOFF = 25;

    public static final double MAX_COPY_NUMBER = 6;
    public static final double CLONAL_COPY_NUMBER_MARGIN = 0.2;

    public static final double DEFAULT_NOISE_READS_PER_MILLION = 30;
    public static final double DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND = 1;

    public static final double LOW_PROBABILITY = 0.05;
    public static final double HIGH_PROBABILITY = 1 - LOW_PROBABILITY;

    public static final int SOMATIC_PEAK_MIN_VARIANTS = 10;
    public static final double SOMATIC_PEAK_MIN_DEPTH_PERC = 0.1;
    public static final int SOMATIC_PEAK_MIN_PEAK_VARIANTS = 5;
    public static final double SOMATIC_PEAK_MIN_PEAK_VARIANTS_PERC = 0.05;

    public static final int SOMATIC_PEAK_MIN_FRAG_VARIANTS = 10;
    public static final int SOMATIC_PEAK_MIN_AVG_DEPTH = 20;
    public static final double SOMATIC_PEAK_NTH_RATIO = 0.08;
    public static final double SOMATIC_PEAK_NTH_RATIO_MIN = 3;
    public static final double SOMATIC_PEAK_BANDWIDTH_MAX = 3;
    public static final double SOMATIC_PEAK_BANDWIDTH_MIN = 0.2;

    public static final int LOW_COUNT_MODEL_MIN_FRAG_VARIANTS = 5;
    public static final int LOW_COUNT_MODEL_MIN_AVG_DEPTH = 50;
    public static final double DROPOUT_RATE_INCREMENT = 0.1;

    public static final double SYNTHETIC_TUMOR_VAF = 0.5;

    public static final int CHIP_MIN_ALLELE_FRAGS = 5;
    public static final double CHIP_MIN_SAMPLE_PERC = 0.33;

    public static final String PURPLE_CTDNA_SOMATIC_VCF_ID = ".purple.somatic.ctdna.";
}
