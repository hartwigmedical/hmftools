package com.hartwig.hmftools.redux;

public class ReduxConstants
{
    public static final String FILE_ID = "redux";

    public static final int DEFAULT_READ_LENGTH = 151;

    // UMIs
    public static final int DEFAULT_MAX_UMI_BASE_DIFF = 1;
    public static final int MAX_IMBALANCED_UMI_BASE_DIFF = 4;
    public static final int MAX_IMBALANCED_UMI_COUNT = 25;
    public static final int MIN_POLYG_UMI_TAIL_LENGTH = 2;
    public static final int MAX_UMI_BASE_DIFF_JITTER_COLLAPSE = 1;

    public static final char DEFAULT_DUPLEX_UMI_DELIM = '_';

    public static final String CONSENSUS_PREFIX = "CNS_";

    // read unmapping
    public static final int UNMAP_MAX_NON_OVERLAPPING_BASES = 10;
    public static final int UNMAP_MIN_SOFT_CLIP = 20;
    public static final int UNMAP_MIN_HIGH_DEPTH = 1000;
    public static final int UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX = 1000;

    public static final int SUPP_ALIGNMENT_SCORE_MIN = 30;

    // SBX-specfic
    public static final double SBX_CONSENSUS_BASE_THRESHOLD = 0.5;

    public static final byte INVALID_BASE_QUAL = -1;

    // consensus building
    public static int CONSENSUS_MAX_DEPTH = 100;

    // base qual recalibration
    public static final int BQR_SAMPLE_SIZE = 2_000_000;
    public static final int BQR_MIN_MAP_QUAL = 50;
    public static final int BQR_CHR_END_BUFFER = 1000000;

    // base quality recalibration
    public static final double BQR_DUAL_AF_LOW = 0.01;
    public static final double BQR_DUAL_AF_HIGH = 0.075;
    public static final int BQR_DUAL_AD = 2;

    public static final double BQR_NON_DUAL_AF_LOW = 0.05;
    public static final double BQR_NON_DUAL_AF_HIGH = 0.125;
    public static final int BQR_NON_DUAL_AD = 3;
}
