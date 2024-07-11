package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.basequal.jitter.JitterModelParams.REPEAT_UNIT_3_PLUS_LABEL;

import java.util.List;

import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;
import com.hartwig.hmftools.sage.filter.SoftFilterConfig;

public class SageConstants
{
    public static final int DEFAULT_MIN_MAP_QUALITY = 10;
    public static final int DEFAULT_MAX_READ_DEPTH = 1000;
    public static final int DEFAULT_MAX_READ_DEPTH_PANEL = 100_000;
    public static final int DEFAULT_SLICE_SIZE = 100_000;
    public static final int REGION_BLOCK_SIZE = 100;
    public static final int DEFAULT_MAX_PARTITION_SLICES = 10;

    public static final int DEFAULT_READ_LENGTH = 151;

    // read context building
    public static final int DEFAULT_FLANK_LENGTH = 10;
    public static final int MIN_CORE_DISTANCE = 2;
    public static final int MAX_REPEAT_LENGTH = 5;
    public static final int MIN_REPEAT_COUNT = 3;
    public static final int OUTER_MIN_REPEAT_COUNT = 6;

    public static final int MIN_SECOND_CANDIDATE_FULL_READS = 3;
    public static final double MIN_SECOND_CANDIDATE_FULL_READS_PERC = 0.25;

    // base quality recalibration
    public static final double BQR_DUAL_AF_LOW = 0.01;
    public static final double BQR_DUAL_AF_HIGH = 0.075;
    public static final int BQR_DUAL_AD = 2;

    public static final double BQR_NON_DUAL_AF_LOW = 0.05;
    public static final double BQR_NON_DUAL_AF_HIGH = 0.125;
    public static final int BQR_NON_DUAL_AD = 3;

    // dual (illumina)/balanced (ultima): sites where [(AF<1% or AD<3) and AF <7.5%]
    // other: sites where [(AF<5% or AD<4) and AF <12.5%]

    public static final int BQR_SAMPLE_SIZE = 2_000_000;
    public static final int DEFAULT_BQR_MIN_MAP_QUAL = 50;

    // read evidence
    public static final int MATCHING_BASE_QUALITY = 20;
    public static final int CORE_LOW_QUAL_MISMATCH_FACTOR = 8;
    public static final int FLANK_LOW_QUAL_MISMATCHES = 3;
    public static final double SC_READ_EVENTS_FACTOR = 12;

    public static final int SC_INSERT_REF_TEST_LENGTH = 12;
    public static final int SC_INSERT_MIN_LENGTH = 5;
    public static final int MIN_INSERT_ALIGNMENT_OVERLAP = 5;

    public static final int MIN_SOFT_CLIP_MIN_BASE_QUAL = 25;
    public static final int MAX_SOFT_CLIP_LOW_QUAL_COUNT = 5;
    public static final double MIN_SOFT_CLIP_HIGH_QUAL_PERC = 0.75;

    public static final int LONG_GERMLINE_INSERT_LENGTH = 11;
    public static final int LONG_GERMLINE_INSERT_READ_VS_REF_DIFF = 2;

    public static final int EVIDENCE_MIN_MAP_QUAL = 1;

    public static final int CHIMERIC_FRAGMENT_LENGTH_MAX = 1000;
    
    // filtering defaults and constants
    public static final int DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY = 0;
    public static final int DEFAULT_HARD_MIN_TUMOR_QUAL = 50;
    public static final double DEFAULT_HARD_MIN_TUMOR_VAF = 0.002;
    public static final int DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT = 2;
    public static final int DEFAULT_FILTERED_MAX_GERMLINE_ALT_SUPPORT = 3;
    public static final double MAX_INDEL_GERMLINE_ALT_SUPPORT = 0.01;

    public static final double HOTSPOT_MIN_TUMOR_VAF_SKIP_QUAL = 0.08;
    public static final int HOTSPOT_MIN_TUMOR_ALT_SUPPORT_SKIP_QUAL = 8;
    public static final int HOTSPOT_MIN_ALT_BASE_QUAL = 150;

    public static final double VAF_PROBABILITY_THRESHOLD = 1e-14;
    public static final double VAF_PROBABILITY_THRESHOLD_HOTSPOT = 1e-9;

    public static final int DEFAULT_MIN_AVG_BASE_QUALITY = 25;
    public static final int DEFAULT_MIN_AVG_BASE_QUALITY_HOTSPOT = 18;

    public static final int MAX_MAP_QUALITY = 60;
    public static final double DEFAULT_MQ_RATIO_FACTOR = 0; // ie disabled,  but for germline should be set to 2.5
    public static final double MQ_RATIO_SMOOTHING = 3;

    public static final double MAX_READ_EDGE_DISTANCE_PERC = 0.33;
    public static final double MAX_READ_EDGE_DISTANCE_PROB = 0.001;
    public static final int MAX_MAP_QUAL_ALT_VS_REF = 15;

    // variant deduplication
    public static final double INDEL_DEDUP_MIN_MATCHED_LPS_PERCENT = 0.1;

    public static final double STRAND_BIAS_CHECK_THRESHOLD = 0.15;
    public static final double STRAND_BIAS_NON_ALT_MIN_DEPTH = 5;
    public static final double STRAND_BIAS_NON_ALT_MIN_BIAS = 0.25;

    public static final int DOUBLE_JITTER_REPEAT_COUNT = 11;
    public static final int MSI_JITTER_MAX_REPEAT_CHANGE = 5;
    public static final double MSI_JITTER_DEFAULT_ERROR_RATE = 0.0001;
    public static final double MSI_JITTER_NOISE_RATE = 0.00025;
    public static final double MSI_JITTER_HARD_FILTER_NOISE_RATE = 0.05;

    public static final int REQUIRED_UNIQUE_FRAG_COORDS = 3;

    // quality calcs
    public static final int DEFAULT_JITTER_MIN_REPEAT_COUNT = 3;
    public static final double JITTER_QUAL_BOOST_MAX_PERC = 1.3;
    public static final int DEFAULT_BASE_QUAL_FIXED_PENALTY = 12;

    public static final int READ_EDGE_PENALTY_0 = 15;
    public static final int READ_EDGE_PENALTY_1 = 5;

    public static final int DEFAULT_MAP_QUAL_FIXED_PENALTY = 0;
    public static final int DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY = 15;
    public static final double DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY = 7;
    public static final double MAP_QUAL_FACTOR_FIXED_PENALTY = 25;
    public static final int MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY = 10;


    // defaults when in high-depth mode
    public static final int DEFAULT_HIGH_DEPTH_BASE_QUAL = 30;
    public static final int DEFAULT_HIGH_DEPTH_MAP_QUAL_FIXED_PENALTY = -15;
    public static final double DEFAULT_HIGH_DEPTH_MAP_QUAL_RATIO_FACTOR = 2.5;

    public static final int VIS_VARIANT_BUFFER = 200;

    public static final SoftFilterConfig DEFAULT_HOTSPOT_FILTER = new SoftFilterConfig(
            "hotspot", 1e-3, 1e-3, 0.01,
            0, 0, 0.1);

    public static final SoftFilterConfig DEFAULT_PANEL_FILTER = new SoftFilterConfig(
            "panel", 1e-5, 1e-5, 0.02,
            0, 0,0.04);

    public static final SoftFilterConfig DEFAULT_HIGH_CONFIDENCE_FILTER = new SoftFilterConfig(
            "high_confidence", 1e-8, 1e-8, 0.025,
            10, 6,0.04);

    public static final SoftFilterConfig DEFAULT_LOW_CONFIDENCE_FILTER = new SoftFilterConfig(
            "low_confidence", 1e-14, 1e-14, 0.025,
            10, 6,
            0.04);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_HP = new JitterModelParams(
            "A/C/G/T", 0.05, 0.09, 0.13, 0.04,
            -0.11, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_DN = new JitterModelParams(
            "AC/AG/AT/CA/CG/CT/GA/GC/GT/TA/TC/TG", 0.11, 0.15, 0.19, 0.04,
            -0.05, 1);

    public static final JitterModelParams DEFAULT_JITTER_PARAMS_3P = new JitterModelParams(
            REPEAT_UNIT_3_PLUS_LABEL, 0.15, 0.19, 0.23, 0.04,
            -0.01, 1);

    public static final List<JitterModelParams> DEFAULT_JITTER_PARAMS = List.of(
            DEFAULT_JITTER_PARAMS_HP, DEFAULT_JITTER_PARAMS_DN, DEFAULT_JITTER_PARAMS_3P);
}
