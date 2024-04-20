package com.hartwig.hmftools.sage.vcf;

public final class VcfTags
{
    public static final String VERSION_META_DATA = "sageVersion";
    public static final String READ_CONTEXT = "RC";
    public static final String READ_CONTEXT_LEFT_FLANK = "RC_LF";
    public static final String READ_CONTEXT_RIGHT_FLANK = "RC_RF";
    public static final String READ_CONTEXT_INDEX = "RC_IDX";

    public static final String DEDUP_MNV_FILTER = "dedupMnv";
    public static final String DEDUP_MIXED_GERMLINE_SOMATIC_FILTER = "dedupMixedGermlineSomatic";
    public static final String DEDUP_SNV_MNV_FILTER = "dedupSnvMnv";
    public static final String DEDUP_INDEL_FILTER = "dedupIndel";
    public static final String DEDUP_MATCH = "dedupMatch";

    public static final String READ_CONTEXT_JITTER = "RC_JIT";
    public static final String READ_CONTEXT_JITTER_DESC = "Read context jitter [Shortened, Lengthened, QualityPenalty]";

    public static final String READ_CONTEXT_EVENTS = "RC_NM";
    public static final String READ_CONTEXT_EVENTS_DESC = "Minimum number of events in read";

    public static final String READ_CONTEXT_REPEAT_SEQUENCE = "RC_REPS";
    public static final String READ_CONTEXT_REPEAT_SEQUENCE_DESC = "Repeat sequence at read context";

    public static final String READ_CONTEXT_MICRO_HOMOLOGY = "RC_MH";
    public static final String READ_CONTEXT_MICRO_HOMOLOGY_DESC = "Micro-homology at read context";

    public static final String READ_CONTEXT_AF_DESC =
            "Allelic frequency calculated from read context counts as (Full + Partial + Core + Realigned + Alt) / Coverage";

    public static final String READ_CONTEXT_IMPROPER_PAIR = "RC_IPC";
    public static final String READ_CONTEXT_IMPROPER_PAIR_DESC = "Read context improper pair count";

    public static final String RAW_DEPTH = "RDP";
    public static final String RAW_SUPPORT_DEPTH = "RAD";
    public static final String RAW_SUPPORT_BASE_QUALITY = "RABQ";

    public static final String AVG_MAP_QUALITY = "AMQ";
    public static final String AVG_MAP_QUALITY_DESC = "Average map quality count (all,alt)";

    public static final String MAX_READ_EDGE_DISTANCE = "MED";
    public static final String MAX_READ_EDGE_DISTANCE_DESC = "Max read edge distance";

    public static final String FRAG_STRAND_BIAS = "SB";
    public static final String FRAG_STRAND_BIAS_DESC = "Fragment strand bias - percentage of forward-orientation fragments (ref,alt)";

    public static final String READ_STRAND_BIAS = "RSB";
    public static final String READ_STRAND_BIAS_DESC = "Read strand bias - percentage of forward-orientation reads (ref,alt)";

    public static final String AVG_BASE_QUAL = "ABQ";
    public static final String AVG_BASE_QUAL_DESC = "Average calculated base quality";

    public static final String MIXED_SOMATIC_GERMLINE = "MSG";
    public static final String MIXED_SOMATIC_GERMLINE_DESC = "Mixed Somatic and Germline variants";

    public static final String LOCAL_PHASE_SET_READ_COUNT = "LPS_RC";
    public static final String LPS_READ_COUNT_DESC = "Local Phase Set Read Count";

    public static final String QUAL_MODEL_TYPE = "QMT";
    public static final String QUAL_MODEL_TYPE_DESC = "Qual model type if applicable";

    public static final String TOTAL_RAW_BASE_QUAL = "TRBQ";
    public static final String TOTAL_RAW_BASE_QUAL_DESC = "Total raw base qual";
}
