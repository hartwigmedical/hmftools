package com.hartwig.hmftools.amber;

public class AmberConstants
{
    public static final String APP_NAME = "Amber";

    public static final int DEFAULT_MIN_BASE_QUALITY = 13;
    public static final int DEFAULT_MIN_MAPPING_QUALITY = 50;
    public static final int DEFAULT_TUMOR_ONLY_MIN_SUPPORT = 2;
    public static final double DEFAULT_TUMOR_ONLY_MIN_VAF = 0.025;

    public static final int DEFAULT_TUMOR_ONLY_MIN_DEPTH = 25;
    public static final int DEFAULT_TUMOR_MIN_DEPTH = 8;

    public static final double DEFAULT_MIN_DEPTH_PERCENTAGE = 0.5;
    public static final double DEFAULT_MAX_DEPTH_PERCENTAGE = 1.5;
    public static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    public static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    public static final int MIN_NORMAL_READ_DEPTH = 7;

    public static final int THREE_PLUS_READS_MIN = 10;
    public static final double THREE_PLUS_READS_SITE_PERC = 0.03;
    public static final double THREE_PLUS_READS_SITE_LOW_VAF_PERC = 0.002;
    public static final double THREE_PLUS_READS_VAF_MIN = 0.05;

    public static final double QUAL_FILTERED_THRESHOLD = 0.15;

    public static final int HOMOZYGOUS_REGION_MIN_SIZE = 500_000;
    public static final int HOMOZYGOUS_REGION_MIN_SNP_LOCI_COUNT = 50;
    public static final int HOMOZYGOUS_REGION_WINDOW_SIZE = 50;
    public static final int HOMOZYGOUS_REGION_MAX_HET_IN_WINDOW = 5;
    public static final double HOMOZYGOUS_REGION_MAX_HET_RATIO = 0.05;
    public static final int HOMOZYGOUS_REGION_LONG_SIZE = 3_000_000;
    public static final int UNIPARENTAL_DISOMY_MIN_LENGTH = 10_000_000;

    public static final int CRAM_MIN_GAP_START = 10000;
    public static final int BAM_MIN_GAP_START = 4000;

    public static final int TARGET_REGION_SITE_BUFFER = 300;

    // Purity analysis
    public static final double LOWER_CDF_BOUND_FOR_CAPTURE = 0.16;
    public static final double UPPER_CDF_BOUND_FOR_CAPTURE = 0.84;
    public static final int MINIMUM_CAPTURED_POINTS = 15;
    public static final double GNOMAD_FREQUENCY_TOLERANCE = 0.15;
    public static final double PEAK_SEARCH_START = 0.005;
    public static final double PEAK_SEARCH_END = 0.37;
    public static final double PEAK_SEARCH_STEP_RATIO = 1.05;
    public static final double PEAK_SEARCH_INITIAL_STEP = 0.001;
    public static final double PEAK_SEARCH_OVERSHOOT = 0.0001;
    public static final double HET_VAF_LOWER_BOUND = 0.35;
    public static final double HET_VAF_UPPER_BOUND = 0.65;
    public static final double MIN_CUTOFF = 0.04;
    public static final double HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION = 0.25;
    public static final double HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION = 0.75;
    public static final double MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION = 0.825;
    public static final double CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION = 0.5;
    public static final double CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS = 0.1;
}
