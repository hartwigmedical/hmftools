package com.hartwig.hmftools.sage.vcf;

public final class VcfTags
{
    public static final String VERSION_META_DATA = "sageVersion";

    // v3.5 onwards
    public static final String READ_CONTEXT_INFO = "RC_INFO";
    public static final String READ_CONTEXT_INFO_DESC = "Read context: alignment start, variant index, left-flank, core, right-flank, read cigar";

    // v3.4 and earlier
    public static final String READ_CONTEXT_CORE = "RC";
    public static final String READ_CONTEXT_LEFT_FLANK = "RC_LF";
    public static final String READ_CONTEXT_RIGHT_FLANK = "RC_RF";
    public static final String READ_CONTEXT_INDEX = "RC_IDX";

    public static final String READ_CONTEXT_UPDATED = "RC_UPDATED";
    public static final String READ_CONTEXT_UPDATED_DESC = "Read context was updated from an earlier Sage version";

    public static final String READ_CONTEXT_JITTER = "RC_JIT";
    public static final String READ_CONTEXT_JITTER_DESC = "Read context jitter [Shortened, Lengthened, QualityPenalty]";

    public static final String READ_CONTEXT_EVENTS = "RC_NM";
    public static final String READ_CONTEXT_EVENTS_DESC = "Minimum number of events in read";

    public static final String READ_CONTEXT_AF_DESC =
            "Allelic frequency calculated from read context counts as (Full + Partial + Core + Realigned + Alt) / Coverage";

    public static final String READ_CONTEXT_IMPROPER_PAIR = "RC_IPC";
    public static final String READ_CONTEXT_IMPROPER_PAIR_DESC = "Read context improper pair count";

    public static final String AVG_MAP_QUALITY = "AMQ";
    public static final String AVG_MAP_QUALITY_DESC = "Average map quality count (all,alt)";

    public static final String MAX_READ_EDGE_DISTANCE = "MED";
    public static final String MAX_READ_EDGE_DISTANCE_DESC = "Max read edge distance";

    public static final String AVG_READ_EDGE_DISTANCE = "AED";
    public static final String AVG_READ_EDGE_DISTANCE_DESC = "Average read edge distance [alt,total]";

    public static final String FRAG_STRAND_BIAS = "SB";
    public static final String FRAG_STRAND_BIAS_DESC = "Fragment strand bias - percentage of forward-orientation fragments (ref,alt)";

    public static final String READ_STRAND_BIAS = "RSB";
    public static final String READ_STRAND_BIAS_DESC = "Read strand bias - percentage of forward-orientation reads (ref,alt)";

    public static final String AVG_BASE_QUAL = "ABQ";
    public static final String AVG_BASE_QUAL_DESC = "Average calculated base quality (all,alt)";

    public static final String AVG_MODIFIED_BASE_QUAL = "AMBQ";
    public static final String AVG_MODIFIED_BASE_QUAL_DESC = "Average modified map quality";

    public static final String AVG_MODIFIED_ALT_MAP_QUAL = "AMMQ";
    public static final String AVG_MODIFIED_ALT_MAP_QUAL_DESC = "Average alt support modified map quality";

    public static final String MIXED_SOMATIC_GERMLINE = "MSG";
    public static final String MIXED_SOMATIC_GERMLINE_DESC = "Mixed Somatic and Germline variants";

    public static final String LOCAL_PHASE_SET_READ_COUNT = "LPS_RC";
    public static final String LPS_READ_COUNT_DESC = "Local Phase Set Read Count";

    public static final String QUAL_MODEL_TYPE = "QMT";
    public static final String QUAL_MODEL_TYPE_DESC = "Qual model type if applicable";

    public static final String TUMOR_QUALITY_PROB = "TQP";
    public static final String TUMOR_QUALITY_PROB_DESC = "Probability as used in min tumor quality filter";

    public static final String MAP_QUAL_FACTOR = "MQF";
    public static final String MAP_QUAL_FACTOR_DESC = "Map qual heuristic as used in min tumor quality filter";

    public static final String SIMPLE_ALT_COUNT = "SAC";
    public static final String SIMPLE_ALT_COUNT_DESC = "Simple alt match count";
}
