package com.hartwig.hmftools.purple.config;

public class PurpleConstants
{
    // constants prefixed with 'DEFAULT' can be overridden in config

    // common
    public static final int WINDOW_SIZE = 1000;

    // no-tumor
    public static final int NO_TUMOR_BAF_TOTAL = 3000;
    public static final double NO_TUMOR_DEPTH_RATIO_MIN = 0.8;
    public static final double NO_TUMOR_DEPTH_RATIO_MAX = 1.2;

    // somatic fitting
    public static final double SNV_HOTSPOT_VAF_PROBABILITY = 0.01;
    public static final int SNV_HOTSPOT_MAX_SNV_COUNT = 2000;

    // SV recovery
    public static final int DEFAULT_RECOVERY_MIN_MATE_QUAL_SCORE = 300;
    public static final int DEFAULT_RECOVERY_MIN_SGL_QUAL_SCORE = 800;
    public static final double RECOVERY_MIN_LENGTH = 1000;
    public static final double RECOVERY_MIN_PLOIDY = 0.5;
    public static final double RECOVERY_MIN_PLOIDY_PERC = 0.5;
    public static final int RECOVERY_MIN_MATE_UNCERTAINTY = 150;
    public static final int RECOVERY_UNBALANCED_MIN_DEPTH_WINDOW_COUNT = 5;
    public static final double RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE = 0.6;
    public static final double RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_PERC = 0.2;

    // germline deletions
    public static final int GERMLINE_DEL_GENE_BUFFER = 500;
    public static final int GERMLINE_DEL_REGION_MIN = 1000;
    public static final int GERMLINE_DEL_REGION_MATCH_BUFFER = 2000;
    public static final double GERMLINE_DEL_CN_CONSISTENCY_MIN = 0.5;
    public static final double GERMLINE_DEL_CN_CONSISTENCY_MACN_PERC = 0.2;
    public static final double GERMLINE_DEL_NORMAL_RATIO = 0.65;
    public static final int GERMLINE_DEL_COHORT_FREQ = 4;
}
