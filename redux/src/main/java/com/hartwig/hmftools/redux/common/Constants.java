package com.hartwig.hmftools.redux.common;

public class Constants
{
    // BAM processing
    public static final int DEFAULT_PARTITION_SIZE = 1000000;
    public static final int DEFAULT_POS_BUFFER_SIZE = 500;

    public static final int DEFAULT_READ_LENGTH = 151;

    // UMIs
    public static final int DEFAULT_MAX_UMI_BASE_DIFF = 1;
    public static final int MAX_IMBALANCED_UMI_BASE_DIFF = 4;
    public static final int MAX_IMBALANCED_UMI_COUNT = 25;

    public static final char DEFAULT_DUPLEX_UMI_DELIM = '_';

    public static final String CONSENSUS_PREFIX = "CNS_";

    // read unmapping
    public static final int UNMAP_MAX_NON_OVERLAPPING_BASES = 10;
    public static final int UNMAP_MIN_SOFT_CLIP = 20;
    public static final int UNMAP_MIN_HIGH_DEPTH = 1000;
    public static final int UNMAP_CHIMERIC_FRAGMENT_LENGTH_MAX = 1000;

    public static final int SUPP_ALIGNMENT_SCORE_MIN = 30;

    // performance
    public static final double LOCK_ACQUIRE_LONG_TIME_MS = 100;

    // consensus building
    public static int CONSENSUS_MAX_DEPTH = 100;
}
