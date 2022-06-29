package com.hartwig.hmftools.svtools.sv_prep;

public final class SvConstants
{
    // region processing
    public static final int DEFAULT_CHR_PARTITION_SIZE = 1000000;
    public static final int DEFAULT_BUCKET_SIZE = 1000;

    // read filtering
    public static final int MIN_ALIGNMENT_BASES = 50;
    public static final int MIN_MAP_QUALITY = 20;
    public static final int MIN_INSERT_ALIGNMENT_OVERLAP = 5;
    public static final int MIN_SOFT_CLIP_LENGTH = 30;
    public static final int MIN_SOFT_CLIP_MIN_BASE_QUAL = 25;
    public static final double MIN_SOFT_CLIP_HIGH_QUAL_PERC = 0.85;
    public static final int MIN_SUPPORTING_READ_DISTANCE = 50;
    public static final int MIN_INSERT_LENGTH_SUPPORT = 10;

    public static final int DEFAULT_FRAG_LENGTH_MIN = 200;
    public static final int DEFAULT_FRAG_LENGTH_MAX = 1500;

    // to confirm
    public static final short DEFAULT_SINGLE_MAP_QUALITY = 255;
    public static short MULTI_MAP_QUALITY_THRESHOLD = 3; // multi-mapped fragments are given map quals of 3 or lower

}
