package com.hartwig.hmftools.cobalt;

public class CobaltConstants
{
    public static final String APP_NAME = "Cobalt";

    public static final int INVALID_VALUE_INDICATOR = -1;

    public static final double DEFAULT_GC_RATIO_MIN = 0.2;
    public static final double DEFAULT_GC_RATIO_MAX = 0.6;

    public static double GC_RATIO_MIN = DEFAULT_GC_RATIO_MIN;
    public static double GC_RATIO_MAX = DEFAULT_GC_RATIO_MAX;

    public static final int WINDOW_SIZE = 1000;
    public static final int PARTITION_SIZE = 100_000_000;
    public static final double MIN_OFF_TARGET_WINDOW_RATIO = 0.5;

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 10;
    public static final int DEFAULT_PCF_GAMMA = 100;

    public static final int ROLLING_MEDIAN_MAX_DISTANCE = 5000;
    public static final int ROLLING_MEDIAN_MIN_COVERAGE = 1000;

    // when we consolidate sparse windows we use this to avoid consolidating over centromere
    public static final int MAX_SPARSE_CONSOLIDATE_DISTANCE = 3_000_000;
}
