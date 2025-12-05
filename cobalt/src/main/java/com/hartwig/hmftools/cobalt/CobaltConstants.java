package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.common.genome.gc.GCBucket.calcGcBucket;

public class CobaltConstants
{
    public static final String APP_NAME = "Cobalt";

    public static final double DEFAULT_GC_RATIO_MIN = 0.24;
    public static final double DEFAULT_GC_RATIO_MAX = 0.68;

    public static double GC_RATIO_MIN = DEFAULT_GC_RATIO_MIN;
    public static double GC_RATIO_MAX = DEFAULT_GC_RATIO_MAX;

    public static int GC_BUCKET_MIN = calcGcBucket(DEFAULT_GC_RATIO_MIN);
    public static int GC_BUCKET_MAX = calcGcBucket(DEFAULT_GC_RATIO_MAX);

    public static final int WINDOW_SIZE = 1000;
    public static final int PARTITION_SIZE = 100_000_000;

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 10;
    public static final int DEFAULT_PCF_GAMMA = 100;

    public static final int ROLLING_MEDIAN_MAX_DISTANCE = 5000;
    public static final int ROLLING_MEDIAN_MIN_COVERAGE = 1000;

    // when we consolidate sparse windows we use this to avoid consolidating over centromere
    public static final int MAX_SPARSE_CONSOLIDATE_DISTANCE = 3_000_000;
}
