package com.hartwig.hmftools.isofox.fusion;

public class FusionConstants
{
    public final static int REALIGN_MIN_SOFT_CLIP_BASE_LENGTH = 3;
    public final static int REALIGN_MAX_SOFT_CLIP_BASE_LENGTH = 10;

    public static final int JUNCTION_BASE_LENGTH = 10; // bases to record from the ref genome around the fusion junction

    public final static int SOFT_CLIP_JUNC_BUFFER = 3; // max that a realigned fragment's position can overhang the fusion junction

}
