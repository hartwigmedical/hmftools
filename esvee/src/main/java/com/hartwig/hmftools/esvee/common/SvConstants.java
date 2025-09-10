package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;

import com.hartwig.hmftools.common.sequencing.SequencingType;

public final class SvConstants
{
    public static final String BAM_HEADER_SAMPLE_INDEX_TAG = "sampleId";

    // commonly used thresholds
    public static final int MIN_VARIANT_LENGTH = 32;
    public static final int MIN_MAP_QUALITY = 20;

    public static final int DEFAULT_MAX_CONCORDANT_FRAG_LENGTH = 1000; // default, otherwise set from BAM fragment sampling
    public static final int MIN_UPPER_FRAGMENT_LENGTH = 800; // in place for panels to maintain a minimum

    public static int maxConcordantFragmentLength(int observedMaxFragmentLength)
    {
        return max(observedMaxFragmentLength, MIN_UPPER_FRAGMENT_LENGTH);
    }

    // sequencing type
    public static SequencingType Sequencing = SequencingType.ILLUMINA;

    // indels
    public static final int MIN_INDEL_SUPPORT_LENGTH = 3;
    public static final int MIN_INDEL_LENGTH = MIN_VARIANT_LENGTH;

    // qual calcs and filters
    public static final double QUAL_CALC_FRAG_SUPPORT_FACTOR = 5;

    public static final int MIN_ANCHOR_LENGTH = 50;

    // LINE elements
    public static final int LINE_MIN_EXTENSION_LENGTH = LINE_POLY_AT_REQ;
    public static final int LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH = LINE_MIN_EXTENSION_LENGTH / 2;
    public static final double LINE_REF_BASE_REPEAT_FACTOR = 1.5;
    public static final int LINE_INDEL_MAX_OVERLAP = 40;
    public static final int LINE_INDEL_MAX_GAP = 30;
}
