package com.hartwig.hmftools.amber;

public class AmberConstants
{
    public static final String APP_NAME = "Amber";

    public static final int DEFAULT_MIN_BASE_QUALITY = 30;
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

}
