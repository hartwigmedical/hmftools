package com.hartwig.hmftools.linx.types;

public class SvaConstants
{
    public static int MIN_TEMPLATED_INSERTION_LENGTH = 30;
    public static int NO_DB_MARKER = -1000; // bigger than longest assembly read distance
    public static int MIN_DEL_LENGTH = 32;

    public static int SHORT_TI_LENGTH = 1000;

    public static int MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE = 5000000;

    public static int MAX_FOLDBACK_CHAIN_LENGTH = 5000;

    public static double LOW_CN_CHANGE_SUPPORT = 0.5;

    public static double SUBCLONAL_LOW_CNC_PERCENT = 0.5;

    public static int MAX_SV_REPLICATION_MULTIPLE = 32;
    public static int MAX_CLUSTER_COUNT_REPLICATION = 500;

    // exclude clusters with too many SVs from chaining
    public static int DEFAULT_CHAINING_SV_LIMIT = 2000;
}
