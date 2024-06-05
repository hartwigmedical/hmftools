package com.hartwig.hmftools.common.variant;

public final class SageVcfTags
{
    public static final String TIER = "TIER";
    public static final String TIER_DESC = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";

    public static final String LOCAL_PHASE_SET = "LPS";
    public static final String LOCAL_PHASE_SET_DESC = "Local Phase Set";

    // NOTE: most downstream applications use reference and read-context repeat and homology information
    public static final String REPEAT_COUNT = "REP_C";
    public static final String REPEAT_COUNT_DESC = "Repeat sequence count";

    public static final String REPEAT_SEQUENCE = "REP_S";
    public static final String REPEAT_SEQUENCE_DESC = "Repeat sequence";

    public static final String READ_CONTEXT_REPEAT_COUNT = "RC_REPC";
    public static final String READ_CONTEXT_REPEAT_COUNT_DESC = "Repeat count from read context";

    public static final String READ_CONTEXT_REPEAT_SEQUENCE = "RC_REPS";
    public static final String READ_CONTEXT_REPEAT_SEQUENCE_DESC = "Repeat sequence from read context";

    public static final String MICROHOMOLOGY = "MH";
    public static final String MICROHOMOLOGY_DESC = "Microhomology";

    public static final String READ_CONTEXT_MICROHOMOLOGY = "RC_MH";
    public static final String READ_CONTEXT_MICROHOMOLOGY_DESC = "Microhomology from read context";

    public static final String TRINUCLEOTIDE_CONTEXT = "TNC";
    public static final String TRINUCLEOTIDE_CONTEXT_DESC = "Tri-nucleotide context";

    public static final String READ_CONTEXT_COUNT = "RC_CNT";
    public static final String READ_CONTEXT_COUNT_DESC =
            "Read context counts [Full, PartialCore, Core, Realigned, Reference, Total]";

    public static final String READ_CONTEXT_QUALITY = "RC_QUAL";
    public static final String READ_CONTEXT_QUALITY_DESC =
            "Read context quality [Full, PartialCore, Core, Realigned, Reference, Total]";

    public static final String UMI_TYPE_COUNTS = "UMI_CNT";
    public static final String UMI_TYPE_COUNTS_DESC =
            "UMI type counts [TotalNone,TotalSingle,TotalDualStrand,AltNone,AltSingle,AltDualStrand]";

    public static final int UMI_TYPE_COUNT = 6;

    public static final String LIST_SEPARATOR = ",";
}
