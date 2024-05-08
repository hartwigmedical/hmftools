package com.hartwig.hmftools.common.sv.gridss;

public class GridssVcfTags
{
    // Gridss-only tags
    public static final String CIRPOS = "CIRPOS";
    public static final String CIRPOS_DESC = "CIRPOS";

    public static final String EVENT = "EVENT";
    public static final String EVENT_DESC = "Linking SV event ID";

    public static final String BEID = "BEID";
    public static final String BEID_DESC = "Identifiers of assemblies supporting the variant";

    public static final String BEIDL = "BEIDL";
    public static final String BEIDL_DESC = "Local alignment offset of corresponding BEID assembly";

    public static final String BEIDH = "BEIDL";
    public static final String BEIDH_DESC = "Remote alignment offset of corresponding BEID assembly";

    public static final String BEALN = "BEALN";
    public static final String BEALN_DESC = "Potential alignment locations of breakend sequence in the format chr:start|strand|cigar|mapq";

    public static final String SPLIT_READS = "SR";
    public static final String SPLIT_READS_DESC = "Count of split reads supporting breakpoint";

    public static final String DISCORDANT_READS = "DP";
    public static final String DISCORDANT_READS_DESC = "Count of discordant reads supporting breakpoint";

    public static final String SV_FRAG_COUNT = "VF";
    public static final String SV_FRAG_COUNT_DESC = "Count of fragments supporting the variant breakpoint";

    public static final String ANCHOR_SUPPORT_CIGAR = "SC";
    public static final String ANCHOR_SUPPORT_CIGAR_DESC = "CIGAR of local anchor, one per assembly";

    public static final String ANCHOR_SUPPORT_CIGAR_LENGTH = "SC_LEN";
    public static final String ANCHOR_SUPPORT_CIGAR_LENGTH_DESC = "Length of local anchor CIGAR, one per assembly";


    public static final String IMPRECISE = "IMPRECISE";
    public static final String IMPRECISE_DESC = "Imprecise structural variation";

    public static final String READ_PAIRS = "RP";
    public static final String READ_PAIRS_DESC = "Count of read pairs supporting breakpoint";

    public static final String INDEL_COUNT = "IC";

    public static final String PAR_ID = "PARID";
    public static final String TAF = "TAF";
    public static final String TAF_DESC = "Tumor allelic frequency (fragment support / total support)";

    public static final String INS_SEQ = "SVINSSEQ";
    public static final String INS_SEQ_DESC = "SVINSSEQ";

    public static final String SGL_FRAG_COUNT = "BVF";
    public static final String SGL_FRAG_COUNT_DESC = "BVF";

    public static final String GRIDSS_BQ = "BQ";
    public static final String GRIDSS_BAQ = "BAQ";
    public static final String GRIDSS_SRQ = "SRQ";
    public static final String GRIDSS_RPQ = "RPQ";
    public static final String GRIDSS_BUMQ = "BUMQ";
    public static final String GRIDSS_BUM = "BUM";
    public static final String GRIDSS_ASRP = "ASRP";
    public static final String GRIDSS_ASSR = "ASSR";
    public static final String GRIDSS_BASRP = "BASRP";
    public static final String GRIDSS_BASSR = "BASSR";
    public static final String GRIDSS_AS = "AS";
    public static final String GRIDSS_CAS = "CAS";
    public static final String GRIDSS_RAS = "RAS";
    public static final String GRIDSS_BSC = "BSC";

    // also set by Gripss when run with Gridss
    public static final String LOCAL_LINKED_BY = "LOCAL_LINKED_BY";
    public static final String LOCAL_LINKED_BY_DESC = "Breakend linking information";

    public static final String REMOTE_LINKED_BY = "REMOTE_LINKED_BY";
    public static final String REMOTE_LINKED_BY_DESC = "Partner breakend linking information";

    public static final String LINKED_BY_DELIM = ",";

    public static final String REALIGN = "REALIGN";
    public static final String REALIGN_DESC = "Variant was realigned";

    // may be merged with SV_TYPE? Gripss uses this to set StructuralVariantType
    public static final String EVENT_TYPE = "EVENTTYPE";
    public static final String EVENT_TYPE_DESC = "Structural variant type";

    public static final String ALT_PATH = "ALT_PATH";
    public static final String ALT_PATH_DESC = "Alternate path";
    public static final String RESCUE_INFO = "RESCUED";
    public static final String RESCUE_INFO_DESC = "Partner breakend rescue";
}
