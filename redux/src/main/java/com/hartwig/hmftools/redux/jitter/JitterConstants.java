package com.hartwig.hmftools.redux.jitter;

import java.util.List;

public class MsiJitterConstants
{
    public static final int MIN_FLANKING_BASE_MATCHES = 5;
    public static final int LOW_BASE_QUAL_FLANKING_BASES = 1;

    public static final int MAX_MICROSAT_UNIT_LENGTH = 5;
    public static final int MIN_MICROSAT_UNIT_COUNT = 4;

    public static final int MIN_ADJACENT_MICROSAT_DISTANCE = 5;

    public static final int MIN_PASSING_SITE_READS = 20;

    public static final double ALT_COUNT_FRACTION_INIT = 0.3;

    public static final double ALT_COUNT_FRACTION_STEP = -0.05;

    public static final double MAX_REJECTED_READ_FRACTION = 0.2;

    public static final double MIN_SITE_READS_BEFORE_OUTLIER_CHECK = 3;
    public static final double DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION = 0.8;

    public static final List<String> SINGLE_BASE_1 = List.of("A", "T");
    public static final List<String> SINGLE_BASE_2 = List.of("C", "G");
    public static final List<String> DUAL_BASE_1 = List.of("AT", "TA");
    public static final List<String> DUAL_BASE_2 = List.of("AC", "CA", "GT", "TG");
    public static final List<String> DUAL_BASE_3 = List.of("AG", "GA", "CT", "TC");
    public static final List<String> DUAL_BASE_4 = List.of("CG", "GC");
}
