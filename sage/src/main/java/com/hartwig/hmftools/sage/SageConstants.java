package com.hartwig.hmftools.sage;

import com.hartwig.hmftools.sage.filter.SoftFilterConfig;

public class SageConstants
{
    public static final int DEFAULT_MIN_MAP_QUALITY = 10;
    public static final int DEFAULT_MAX_READ_DEPTH = 1000;
    public static final int DEFAULT_MAX_READ_DEPTH_PANEL = 100_000;
    public static final int DEFAULT_SLICE_SIZE = 100_000;

    public static final int DEFAULT_READ_CONTEXT_FLANK_SIZE = 10;
    public static final int MIN_CORE_DISTANCE = 2;

    public static final int MIN_SECOND_CANDIDATE_FULL_READS = 3;
    public static final double MIN_SECOND_CANDIDATE_FULL_READS_PERC = 0.25;

    public static final boolean DEFAULT_MNV = true;

    // base quality recalibration
    public static final int DEFAULT_BQR_MAX_ALT_COUNT = 3;
    public static final int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;
    public static final int DEFAULT_BQR_MIN_MAP_QUAL = 10;

    public static final int MATCHING_BASE_QUALITY = 20;
    public static final int CORE_LOW_QUAL_MISMATCH_BASE_LENGTH = 20;
    public static final double SC_READ_EVENTS_FACTOR = 12;

    public static final int SC_INSERT_MIN_FLANK_LENGTH = 10;
    public static final int SC_INSERT_MIN_LENGTH = 5;

    public static final int NORMAL_RAW_ALT_BQ_MAX = 25;
    public static final int LONG_GERMLINE_INSERT_LENGTH = 10;

    public static final int DEFAULT_EVIDENCE_MAP_QUAL = 1;
    
    // filtering defaults and constants
    public static final int DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY = 0;
    public static final int DEFAULT_HARD_MIN_TUMOR_QUAL = 50;
    public static final double DEFAULT_HARD_MIN_TUMOR_VAF = 0.01;
    public static final int DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT = 2;
    public static final int DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT = 3;

    public static final double HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL = 0.08;
    public static final int HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL = 8;
    public static final int HOTSPOT_MIN_RAW_ALT_BASE_QUAL = 150;

    public static final SoftFilterConfig DEFAULT_HOTSPOT_FILTER = new SoftFilterConfig(70, 0.005,
            0, 0, 0, 0,
            0.1, 0.5);

    public static final SoftFilterConfig DEFAULT_PANEL_FILTER = new SoftFilterConfig(100, 0.02  ,
            0, 0, 0, 0,
            0.04, 0.04);

    public static final SoftFilterConfig DEFAULT_HIGH_CONFIDENCE_FILTER = new SoftFilterConfig(160, 0.025,
            10, 15, 6, 9,
            0.04, 0.04);

    public static final SoftFilterConfig DEFAULT_LOW_CONFIDENCE_FILTER = new SoftFilterConfig(240, 0.025,
            10, 15, 6, 9,
            0.04, 0.04);


    public static final String ITEM_DELIM = ";";
    public static final String SUB_ITEM_DELIM = ":";
}
