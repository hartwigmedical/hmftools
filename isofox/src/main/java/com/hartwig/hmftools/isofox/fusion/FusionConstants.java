package com.hartwig.hmftools.isofox.fusion;

public class FusionConstants
{
    public static final int REALIGN_MIN_SOFT_CLIP_BASE_LENGTH = 3;
    public static final int REALIGN_MAX_SOFT_CLIP_BASE_LENGTH = 10;

    public static final int JUNCTION_BASE_LENGTH = 10; // bases to record from the ref genome around the fusion junction

    public static final int SOFT_CLIP_JUNC_BUFFER = 3; // max that a realigned fragment's position can overhang the fusion junction

    public static final int HIGH_LOG_COUNT = 10000;
}
