package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.common.genome.gc.GCProfile.MIN_MAPPABLE_PERCENTAGE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.WINDOW_SIZE;

public class NormConstants
{
    public static final int REGION_SIZE = WINDOW_SIZE;

    public static final double MAPPABILITY_THRESHOLD = MIN_MAPPABLE_PERCENTAGE;

    public static final double MIN_ENRICHMENT_RATIO = 0.1;

    public static final double MIN_ADJACENT_ENRICHMENT_RATIO = 0.15;

    public static final double MAX_WGS_ADJUST_RATIO = 0.3; // vs the median
    public static final double MEDIAN_PERCENTILE = 0.5;
    public static final double WGS_PERCENTILE_DAMPEN = 0.2;

    // region filtering
    public static final double FILTER_CHECK_PERCENTILE = 0.2;
    public static final double FILTER_MEDIAN_MAX_RATIO = 3;
    public static final double FILTER_CHECK_PERCENTILE_RATIO_THRESHOLD = 0.2;
    public static final double FILTER_CHECK_PERCENTILE_STD_DEV_THRESHOLD = 0.5;
    public static final int FILTER_COHORT_MIN = 10;
    public static final double FILTER_COHORT_MIN_PERCENT = 0.5;

}
