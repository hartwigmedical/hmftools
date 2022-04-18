package com.hartwig.hmftools.cobalt;

public class CobaltConstants
{
    public static final double MIN_DIPLOID = 0.85;
    public static final double MAX_DIPLOID = 1.15;
    public static final int WINDOW_SIZE = 1000;
    public static final int OFF_TARGET_WINDOW_SIZE = 100_000;
    public static final double MIN_OFF_TARGET_WINDOW_RATIO = 0.5;

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 10;
    public static final int DEFAULT_PCF_GAMMA = 100;

    public static final int ROLLING_MEDIAN_MAX_DISTANCE = 5000;
    public static final int ROLLING_MEDIAN_MIN_COVERAGE = 1000;
}
