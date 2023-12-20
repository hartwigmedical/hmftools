package com.hartwig.hmftools.sage;

import com.hartwig.hmftools.sage.filter.SoftFilterConfig;

public class SageConstants
{
    public static final int DEFAULT_MIN_MAP_QUALITY = 10;
    public static final int DEFAULT_MAX_READ_DEPTH = 1000;
    public static final int DEFAULT_MAX_READ_DEPTH_PANEL = 100_000;
    public static final int DEFAULT_SLICE_SIZE = 100_000;
    public static final int DEFAULT_MAX_PARTITION_SLICES = 10;

    public static final int DEFAULT_READ_LENGTH = 151;
    public static final int DEFAULT_READ_CONTEXT_FLANK_SIZE = 10;
    public static final int MIN_CORE_DISTANCE = 2;

    public static final int MIN_SECOND_CANDIDATE_FULL_READS = 3;
    public static final double MIN_SECOND_CANDIDATE_FULL_READS_PERC = 0.25;

    // base quality recalibration
    public static final double DEFAULT_BQR_MAX_ALT_PERC = 0.05;
    public static final int DEFAULT_BQR_MAX_ALT_COUNT = 3;
    public static final int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;
    public static final int DEFAULT_BQR_MIN_MAP_QUAL = 10;

    // read evidence
    public static final int MATCHING_BASE_QUALITY = 20;
    public static final int CORE_LOW_QUAL_MISMATCH_BASE_LENGTH = 20;
    public static final double SC_READ_EVENTS_FACTOR = 12;
    public static final int REALIGN_READ_MIN_INDEL_LENGTH = 3;
    public static final int REALIGN_READ_CONTEXT_MIN_SEARCH_LENGTH = 20;
    public static final int REALIGN_READ_CONTEXT_MIN_SEARCH_BUFFER = 5;

    public static final int SC_INSERT_MIN_SC_LENGTH = 12;
    public static final int SC_INSERT_MIN_LENGTH = 5;
    public static final int MIN_INSERT_ALIGNMENT_OVERLAP = 5;

    public static final int MIN_SOFT_CLIP_MIN_BASE_QUAL = 25;
    public static final int MAX_SOFT_CLIP_LOW_QUAL_COUNT = 5;
    public static final double MIN_SOFT_CLIP_HIGH_QUAL_PERC = 0.75;

    public static final int NORMAL_RAW_ALT_BQ_MAX = 25;
    public static final int LONG_GERMLINE_INSERT_LENGTH = 10;

    public static final int EVIDENCE_MIN_MAP_QUAL = 1;

    public static final int CHIMERIC_FRAGMENT_LENGTH_MAX = 1000;
    
    // filtering defaults and constants
    public static final int DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY = 0;
    public static final int DEFAULT_HARD_MIN_TUMOR_QUAL = 50;
    public static final double DEFAULT_HARD_MIN_TUMOR_VAF = 0.01;
    public static final int DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT = 2;
    public static final int DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT = 3;
    public static final double MAX_INDEL_GERMLINE_ALT_SUPPORT = 0.01;

    public static final double HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL = 0.08;
    public static final int HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL = 8;
    public static final int HOTSPOT_MIN_RAW_ALT_BASE_QUAL = 150;

    public static final double VAF_PROBABILITY_THRESHOLD = 10e-14;
    public static final double VAF_PROBABILITY_THRESHOLD_HOTSPOT = 10e-9;

    public static final int DEFAULT_MIN_AVG_BASE_QUALITY = 25;
    public static final int DEFAULT_MIN_AVG_BASE_QUALITY_HOTSPOT = 18;

    public static final int MAX_MAP_QUALITY = 60;
    public static final double DEFAULT_MQ_RATIO_FACTOR = 0; // ie disabled,  but for germline should be set to 2.5
    public static final double MQ_RATIO_SMOOTHING = 3;

    public static final int MAX_READ_EDGE_DISTANCE = 40;
    public static final double MAX_READ_EDGE_DISTANCE_PROB = 0.001;

    // variant deduplication
    public static final double INDEL_DEDUP_MIN_MATCHED_LPS_PERCENT = 0.1;

    public static final double STRAND_BIAS_CHECK_THRESHOLD = 0.1;
    public static final double STRAND_BIAS_REF_MIN_DEPTH = 5;
    public static final double STRAND_BIAS_REF_MIN_BIAS = 0.25;

    public static final int JITTER_INDEL_MAX_REPEATS = 3;
    public static final double JITTER_INDEL_VAF_THRESHOLD = 0.015;

    public static final int JITTER_NON_INDEL_MAX_REPEATS = 5;
    public static final double JITTER_NON_INDEL_VAF_THRESHOLD = 0.01;

    public static final int REQUIRED_UNIQUE_FRAG_COORDS = 3;

    // quality calcs
    public static final double DEFAULT_JITTER_PENALTY = 0.25;
    public static final int DEFAULT_JITTER_MIN_REPEAT_COUNT = 3;
    public static final int DEFAULT_BASE_QUAL_FIXED_PENALTY = 12;
    public static final int DEFAULT_READ_EDGE_FIXED_PENALTY = 0;
    public static final int DEFAULT_READ_EDGE_FACTOR = 3;
    public static final int DEFAULT_MAP_QUAL_FIXED_PENALTY = 15;
    public static final int DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY = 15;
    public static final double DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY = 7;
    public static final int DEFAULT_HIGH_DEPTH_BASE_QUAL = 30;

    public static final SoftFilterConfig DEFAULT_HOTSPOT_FILTER = new SoftFilterConfig(
            "hotspot", 55, 0.01,
            0, 0, 0, 0,
            0.1, 0.5);

    public static final SoftFilterConfig DEFAULT_PANEL_FILTER = new SoftFilterConfig(
            "panel", 100, 0.02,
            0, 0, 0, 0,
            0.04, 0.04);

    public static final SoftFilterConfig DEFAULT_HIGH_CONFIDENCE_FILTER = new SoftFilterConfig(
            "high_confidence", 160, 0.025,
            10, 15, 6, 9,
            0.04, 0.04);

    public static final SoftFilterConfig DEFAULT_LOW_CONFIDENCE_FILTER = new SoftFilterConfig(
            "low_confidence", 240, 0.025,
            10, 15, 6, 9,
            0.04, 0.04);
}
