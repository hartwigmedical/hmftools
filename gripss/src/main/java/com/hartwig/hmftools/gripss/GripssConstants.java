package com.hartwig.hmftools.gripss;

public final class GripssConstants
{
    // alternate paths
    public static final int MAX_ASSEMBLY_JUMPS = 5;
    public static final int MAX_TRANSITIVE_JUMPS = 2;
    public static final int MAX_VARIANTS = 500000;
    public static final int MAX_ALTERNATIVES = 25;
    public static final int MAX_ALTERNATIVES_SEEK_DISTANCE = 1000;
    public static final int MAX_ALTERNATIVES_ADDITIONAL_DISTANCE = MAX_ALTERNATIVES_SEEK_DISTANCE;

    public static final int MIN_TRANSITIVE_DISTANCE = 30;
    public static final int MAX_TRANSITIVE_SEEK_DISTANCE = 2000;
    public static final int MAX_TRANSITIVE_ADDITIONAL_DISTANCE = 1000;

    // double-stranded breaks
    public static final int MAX_DSB_SEEK_DISTANCE = 1000;
    public static final int MAX_DSB_DISTANCE = 30;

    // rescue
    public static final int SHORT_RESCUE_LENGTH = 10000;
}
