package com.hartwig.hmftools.markdups.common;

public class Constants
{
    // BAM processing
    public static final int DEFAULT_PARTITION_SIZE = 1000000;
    public static final int DEFAULT_POS_BUFFER_SIZE = 500;

    public static final int DEFAULT_READ_LENGTH = 151;

    // UMIs
    public static final int DEFAULT_MAX_UMI_BASE_DIFF = 1;

    public static final char DEFAULT_DUPLEX_UMI_DELIM = '_';

    // performance
    public static final double LOCK_ACQUIRE_LONG_TIME_MS = 100;
    public static final int BAM_READ_CACHE_BUFFER = 100000;
}
