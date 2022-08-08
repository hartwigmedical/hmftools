package com.hartwig.hmftools.svprep;

import com.hartwig.hmftools.common.sv.LineElements;

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
    public static final int MIN_SOFT_CLIP_LENGTH = 30;
    public static final int MIN_LINE_SOFT_CLIP_LENGTH = LineElements.LINE_POLY_AT_TEST_LEN;
    public static final int MIN_SOFT_CLIP_MIN_BASE_QUAL = 25;
    public static final double MIN_SOFT_CLIP_HIGH_QUAL_PERC = 0.85;
    public static final int MAX_SOFT_CLIP_LOW_QUAL_COUNT = 5;
    public static final int MIN_INSERT_LENGTH_SUPPORT = 10;
    public static final int MIN_INDEL_LENGTH = 32;
    public static final int MIN_INDEL_SUPPORT_LENGTH = 20;
    public static final int LOW_BASE_QUALITY = 20;

    public static final int REPEAT_BREAK_SC_CHECK_LENGTH = 6;
    public static final int REPEAT_BREAK_MATCH_CHECK_LENGTH = 9;
    public static final int REPEAT_BREAK_MIN_MAP_QUAL = 40;
    public static final int REPEAT_BREAK_MIN_SC_LENGTH = 50;

    // supporting reads
    public static final int JUNCTION_SUPPORT_CAP = 0; // no limit
    public static final int MIN_SUPPORTING_READ_DISTANCE = 50;
    public static final int MAX_DISCORDANT_READ_DISTANCE = 1000;
    public static final int MAX_HIGH_QUAL_BASE_MISMATCHES = 1;

    // final junction filtering
    public static final int MIN_HOTSPOT_JUNCTION_SUPPORT = 2;
    public static final int MIN_JUNCTION_SUPPORT = 3;

    // fragment length distribution and filtering
    public static final int FRAG_LENGTH_DIST_SAMPLE_SIZE = 10000;
    public static final int MAX_FRAGMENT_LENGTH = 1100;
    public static final int DEFAULT_READ_LENGTH = 151;
}
