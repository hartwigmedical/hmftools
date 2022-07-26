package com.hartwig.hmftools.svprep;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public final class SvConstants
{
    // region processing
    public static final int DEFAULT_CHR_PARTITION_SIZE = 1000000;
    public static final int DOWN_SAMPLE_FRACTION = 20; // split partition into segments
    public static final int DOWN_SAMPLE_THRESHOLD = 1_500_000; // per partition

    // candidate junction fragments
    public static final int MIN_ALIGNMENT_BASES = 50;
    public static final int MIN_MAP_QUALITY = 20;
    public static final int MIN_INSERT_ALIGNMENT_OVERLAP = 5;
    public static final int MIN_SOFT_CLIP_LENGTH = 20;
    public static final int MIN_SOFT_CLIP_MIN_BASE_QUAL = 25;
    public static final double MIN_SOFT_CLIP_HIGH_QUAL_PERC = 0.85;
    public static final int MIN_INSERT_LENGTH_SUPPORT = 10;
    public static final int MIN_INDEL_LENGTH = 20;
    public static final int LOW_BASE_QUALITY = 20;

    public static final int REPEAT_BREAK_SC_CHECK_LENGTH = 6;
    public static final int REPEAT_BREAK_MATCH_CHECK_LENGTH = 9;
    public static final int REPEAT_BREAK_MIN_MAP_QUAL = 40;
    public static final int REPEAT_BREAK_MIN_SC_LENGTH = 50;

    // supporting reads
    public static final int JUNCTION_SUPPORT_CAP = 0; // no limit
    public static final int MIN_SUPPORTING_READ_DISTANCE = 50;
    public static final int MAX_DISCORDANT_READ_DISTANCE = 1000;

    // final junction filtering
    public static final int MIN_HOTSPOT_JUNCTION_SUPPORT = 2;
    public static final int MIN_JUNCTION_SUPPORT = 3;

    // fragment length distribution and filtering
    public static final int FRAG_LENGTH_DIST_SAMPLE_SIZE = 10000;
    public static final int MAX_FRAGMENT_LENGTH = 1500;
    public static final int DEFAULT_READ_LENGTH = 151;

    public static final ChrBaseRegion EXCLUDED_REGION_1_REF_37 = new ChrBaseRegion("2", 33141260, 33141700);
    public static final ChrBaseRegion EXCLUDED_REGION_1_REF_38 = new ChrBaseRegion("chr2", 32916190, 32916630);
}
