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
    public static final int SNV_HOTSPOT_MAX_SNV_COUNT = 1000;

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



}
