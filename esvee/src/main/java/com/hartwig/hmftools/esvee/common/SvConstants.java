package com.hartwig.hmftools.esvee.common;

public final class SvConstants
{
    public static final String BAM_HEADER_SAMPLE_INDEX_TAG = "sampleId";

    // commonly used thresholds
    public static final int MIN_VARIANT_LENGTH = 32;
    public static final int DEFAULT_DISCORDANT_FRAGMENT_LENGTH = 1000; // default, otherwise set from BAM fragment sampling
    public static int LOW_BASE_QUAL_THRESHOLD = 26;

    // indels
    public static final int MIN_INDEL_SUPPORT_LENGTH = 5;
    public static final int MIN_INDEL_LENGTH = MIN_VARIANT_LENGTH;

    // qual calcs and filters
    public static final double QUAL_CALC_FRAG_SUPPORT_FACTOR = 5;




}
