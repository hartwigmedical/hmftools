package com.hartwig.hmftools.common.basequal.jitter;

public class JitterAnalyserConstants
{
    public static final int MIN_FLANKING_BASE_MATCHES = 5;

    public static final int MAX_MICROSAT_UNIT_LENGTH = 5;
    public static final int MIN_MICROSAT_UNIT_COUNT = 4;

    public static final int MIN_ADJACENT_MICROSAT_DISTANCE = 5;

    public static final double ALT_COUNT_FRACTION_INIT = 0.3;

    public static final double ALT_COUNT_FRACTION_STEP = -0.05;

    public static final double MAX_REJECTED_READ_FRACTION = 0.2;

    public static final double MIN_SITE_READS_BEFORE_OUTLIER_CHECK = 3;
    public static final double MAX_SINGLE_SITE_ALT_CONTRIBUTION = 0.8;
}
