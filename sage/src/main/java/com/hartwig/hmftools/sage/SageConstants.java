package com.hartwig.hmftools.sage;

public class SageConstants
{
    public static final int DEFAULT_THREADS = 2;
    public static final int DEFAULT_MIN_MAP_QUALITY = 10;
    public static final int DEFAULT_MAX_READ_DEPTH = 1000;
    public static final int DEFAULT_MAX_READ_DEPTH_PANEL = 100_000;
    public static final int DEFAULT_MAX_REALIGNMENT_DEPTH = 1000;
    public static final int DEFAULT_SLICE_SIZE = 100_000;
    public static final int DEFAULT_READ_CONTEXT_FLANK_SIZE = 10;
    public static final boolean DEFAULT_MNV = true;

    // base quality recalibration
    public static final int DEFAULT_BQR_MAX_ALT_COUNT = 3;
    public static final int DEFAULT_BQR_SAMPLE_SIZE = 2_000_000;
    public static final int DEFAULT_BQR_MIN_MAP_QUAL = 10;

    public static final int MATCHING_BASE_QUALITY = 25;

    public static final String ITEM_DELIM = ";";
    public static final String SUB_ITEM_DELIM = ":";
}
