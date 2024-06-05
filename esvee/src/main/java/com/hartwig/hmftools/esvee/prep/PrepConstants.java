package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;

import com.hartwig.hmftools.common.sv.LineElements;

public final class PrepConstants
{
    public static final String PREP_JUNCTIONS_FILE_ID = "junctions" + TSV_EXTENSION;
    public static final String PREP_FRAG_LENGTH_FILE_ID = "fragment_lengths" + TSV_EXTENSION;

    // common fields
    public static final String FLD_JUNCTION_FRAGS = "JunctionFrags";
    public static final String FLD_INDEL_JUNCTION = "Indel";
    public static final String FLD_EXACT_SUPPORT_FRAGS = "ExactSupportFrags";
    public static final String FLD_OTHER_SUPPORT_FRAGS = "OtherSupportFrags";
    public static final String FLD_HOTSPOT_JUNCTION = "Hotspot";

    // region processing
    public static final int DEFAULT_CHR_PARTITION_SIZE = 1000000;

    public static final int DEFAULT_READ_LENGTH = 151;

    // candidate junction fragments
    public static final int MIN_ALIGNMENT_BASES = 50;
    public static final int MIN_MAP_QUALITY = 20;
    public static final int MIN_INSERT_ALIGNMENT_OVERLAP = 5;
    public static final int MIN_SOFT_CLIP_LENGTH = 30;
    public static final int MIN_LINE_SOFT_CLIP_LENGTH = LineElements.LINE_POLY_AT_TEST_LEN;
    public static final double MIN_SOFT_CLIP_HIGH_QUAL_PERC = 0.75;
    public static final int MAX_SOFT_CLIP_LOW_QUAL_COUNT = 5;
    public static final int MIN_INSERT_LENGTH_SUPPORT = 10;

    public static final int REPEAT_BREAK_CHECK_LENGTH = 9;
    public static final int REPEAT_BREAK_MIN_MAP_QUAL = 40;
    public static final int REPEAT_BREAK_MIN_SC_LENGTH = 50;

    // supporting reads
    public static final int MIN_SUPPORTING_READ_DISTANCE = 50;
    public static final int UNPAIRED_READ_JUNCTION_DISTANCE = 5;
    public static final int MAX_SUPPORT_FRAGMENT_DISTANCE = 1000;
    public static final int MAX_HIGH_QUAL_BASE_MISMATCHES = 1;
    public static final double MIN_EXACT_BASE_PERC = 0.25;

    // discordant groups
    public static final int DISCORDANT_GROUP_MIN_FRAGMENTS = 3;
    public static final int DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT = 5;
    public static final int DISCORDANT_GROUP_MAX_DISTANCE = 500;

    // final junction filtering
    public static final int MIN_HOTSPOT_JUNCTION_SUPPORT = 1;
    public static final int MIN_JUNCTION_SUPPORT = 2;

    // fragment length distribution and filtering
    public static final int FRAG_LENGTH_DIST_SAMPLE_SIZE = 100000;
    public static final int FRAG_LENGTH_DIST_MIN_QUAL = 60;
    public static final int DISCORDANT_FRAGMENT_LENGTH_MIN = 1100;
    public static final int FRAG_LENGTH_DIST_MAX_LENGTH = 1500;
    public static final double FRAG_LENGTH_DIST_PERCENTILE = 0.9975;
    public static final double FRAG_LENGTH_1_STD_DEV_PERCENTILE = 0.16;

    public static final int DEFAULT_MAX_FRAGMENT_LENGTH = DISCORDANT_FRAGMENT_LENGTH_MIN;

    public static final String BAM_RECORD_SAMPLE_ID_TAG = "SI";
}
