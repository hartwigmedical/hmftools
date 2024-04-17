package com.hartwig.hmftools.esvee.common;

public final class SvConstants
{
    public static final String ESVEE_FILE_ID = "esvee";
    public static final String PREP_FILE_ID = "esvee.prep";

    public static final String FILE_NAME_DELIM = ".";

    public static final String BAM_HEADER_SAMPLE_INDEX_TAG = "sampleId";

    // commonly used thresholds
    public static final int MIN_VARIANT_LENGTH = 32;
    public static final int DEFAULT_DISCORDANT_FRAGMENT_LENGTH = 1000; // default, otherwise set from BAM fragment sampling
    public static int LOW_BASE_QUAL_THRESHOLD = 26;

    // read adjustments
    public static final int POLY_G_TRIM_LENGTH = 4;

    // read filtering
    public static final double LOW_BASE_TRIM_PERC = 0.3;

    // indels
    public static final int MIN_INDEL_SUPPORT_LENGTH = 5;
    public static final int MIN_INDEL_LENGTH = MIN_VARIANT_LENGTH;


}
