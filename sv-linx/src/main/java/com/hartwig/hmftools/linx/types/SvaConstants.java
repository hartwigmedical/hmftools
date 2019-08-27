package com.hartwig.hmftools.linx.types;

public class SvaConstants
{
    public static final int DEFAULT_PROXIMITY_DISTANCE = 5000;

    public static final int MIN_TEMPLATED_INSERTION_LENGTH = 30;
    public static final int NO_DB_MARKER = -1000; // bigger than longest assembly read distance
    public static final int MIN_DEL_LENGTH = 32;

    public static final int SHORT_TI_LENGTH = 1000;

    public static final long MAX_MERGE_DISTANCE = 5000000;

    public static final int MAX_FOLDBACK_CHAIN_LENGTH = 5000;

    public static final int MAX_SIMPLE_DUP_DEL_CUTOFF = 5000000;
    public static final int MIN_SIMPLE_DUP_DEL_CUTOFF = 100000;

    public static final double LOW_CN_CHANGE_SUPPORT = 0.5;

    public static final double LOW_PLOIDY_THRESHOLD = 0.75;
    public static final double PLOIDY_DIFF_THRESHOLD = 0.5;

    public static final double SUBCLONAL_LOW_CNC_PERCENT = 0.5;

    // exclude clusters with too many SVs from chaining
    public static int DEFAULT_CHAINING_SV_LIMIT = 2000;
}
