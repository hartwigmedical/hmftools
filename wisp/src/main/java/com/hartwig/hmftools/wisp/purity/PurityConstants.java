package com.hartwig.hmftools.wisp.purity;

import java.util.List;

import com.google.common.collect.Lists;

public class PurityConstants
{
    public static final double MAX_SUBCLONAL_LIKELIHOOD = 0.5;
    public static final double SUBCLONAL_VCN_THRESHOLD = 0.7;
    public static final double MAX_REPEAT_COUNT = 3;

    public static final int MIN_QUAL_PER_AD = 18;

    public static final double COPY_NUMBER_MAX = 6;
    public static final double COPY_NUMBER_CLONAL_MARGIN = 0.2;
    public static final double COPY_NUMBER_LOD_FACTOR = 0.004;
    public static final double COPY_NUMBER_LOD_CLONAL_FACTOR = 3;

    public static final double DEFAULT_NOISE_READS_PER_MILLION = 30;
    public static final double DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND = 1;

    public static final double LOW_PROBABILITY = 0.01;
    public static final double HIGH_PROBABILITY = 1 - LOW_PROBABILITY;
    public static final int LOD_MAX_ITERATIONS = 40;
    public static final double LOD_MIN_PROB_DIFF_PERC = 0.001;

    // somatic variants
    public static final int SOMATIC_PEAK_MIN_VARIANTS = 10;
    public static final double SOMATIC_PEAK_MIN_DEPTH_PERC = 0.1;
    public static final int SOMATIC_PEAK_MIN_PEAK_VARIANTS = 3;
    public static final double SOMATIC_PEAK_MIN_PEAK_VARIANTS_PERC = 0.03;

    public static final int SOMATIC_PEAK_MIN_FRAG_VARIANTS = 10;
    public static final int SOMATIC_PEAK_MIN_AVG_DEPTH = 20;
    public static final double SOMATIC_PEAK_NTH_RATIO = 0.08;
    public static final double SOMATIC_PEAK_NTH_RATIO_MIN = 3;
    public static final double SOMATIC_PEAK_BANDWIDTH_MAX = 3;
    public static final double SOMATIC_PEAK_BANDWIDTH_MIN = 0.2;
    public static final double SOMATIC_PEAK_BANDWIDTH_ABS_MIN = 0.1;
    public static final double SOMATIC_PEAK_MAX_IMPLIED_TF = 2;

    public static final int LOW_COUNT_MODEL_MIN_FRAG_VARIANTS = 4;
    public static final int LOW_COUNT_MODEL_MIN_AVG_DEPTH = 50;

    public static final int LOW_COUNT_MODEL_MIN_2_PLUS_FRAGS = 2;
    public static final double LOW_COUNT_MODEL_MIN_2_PLUS_FRAG_PERC = 0.002;

    public static final double DROPOUT_RATE_INCREMENT = 0.1;

    public static final List<Integer> SNV_QUAL_THRESHOLDS = Lists.newArrayList(0, 38, 42, 45);
    public static final int DEFAULT_BQR_MIN_QUAL = 36;
    public static final double BQR_MIN_ERROR_RATE = 1e-5;

    public static final double SYNTHETIC_TUMOR_VAF = 0.5;

    public static final int CHIP_MIN_ALLELE_FRAGS = 5;
    public static final double CHIP_MIN_SAMPLE_PERC = 0.33;

    public static final double AMBER_LOH_MINOR_ALLELE_THRESHOLD = 0.2;
    public static final double AMBER_LOH_CN_THRESHOLD = 0.8;
    public static final double AMBER_LOH_MIN_TUMOR_BAF = 0.55;
    public static final double AMBER_LOH_MIN_AF = 0.55;

    public static final String PURPLE_APPENDED_SOMATIC_VCF_ID = ".purple.somatic.ctdna.";
}
