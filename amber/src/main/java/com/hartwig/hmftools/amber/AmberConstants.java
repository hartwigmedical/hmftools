package com.hartwig.hmftools.amber;

public class AmberConstants
{
    public static final String APP_NAME = "Amber";

    public static final int DEFAULT_MIN_BASE_QUALITY = 13;
    public static final int DEFAULT_MIN_PARTITION = 10000;
    public static final int DEFAULT_MIN_MAPPING_QUALITY = 50;
    public static final int DEFAULT_TUMOR_ONLY_MIN_SUPPORT = 2;
    public static final double DEFAULT_TUMOR_ONLY_MIN_VAF = 0.05;

    public static final int DEFAULT_TUMOR_ONLY_MIN_DEPTH = 25;
    public static final int DEFAULT_TUMOR_MIN_DEPTH = 8;

    public static final double DEFAULT_MIN_DEPTH_PERCENTAGE = 0.5;
    public static final double DEFAULT_MAX_DEPTH_PERCENTAGE = 1.5;
    public static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    public static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    public static final int MIN_NORMAL_READ_DEPTH = 7;
    public static final long MIN_THREE_PLUS_READS = 2000;

    public static final int HOMOZYGOUS_REGION_MIN_SIZE = 500_000;
    public static final int HOMOZYGOUS_REGION_MIN_SNP_LOCI_COUNT = 50;
    public static final int HOMOZYGOUS_REGION_WINDOW_SIZE = 50;
    public static final int HOMOZYGOUS_REGION_MAX_HET_IN_WINDOW = 5;
    public static final double HOMOZYGOUS_REGION_MAX_HET_RATIO = 0.05;
    public static final int HOMOZYGOUS_REGION_LONG_SIZE = 3_000_000;
    public static final int UNIPARENTAL_DISOMY_MIN_LENGTH = 10_000_000;

    // BAM/CRAM slicing optimisations
    public static final int BAM_REGION_GROUP_MAX = 150000;
    public static final int CRAM_REGION_GROUP_MAX = 10000;

    public static final int CRAM_MIN_GAP_START = 10000;
    public static final int BAM_MIN_GAP_START = 4000;

    public static final int CRAM_MIN_GAP_INCREMENT = 1000;
    public static final int BAM_MIN_GAP_INCREMENT = 200;

    public static final int TARGET_REGION_SITE_BUFFER = 300;

}
