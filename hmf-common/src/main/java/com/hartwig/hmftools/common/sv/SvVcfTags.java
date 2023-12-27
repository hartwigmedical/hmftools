package com.hartwig.hmftools.common.sv;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class SvVcfTags
{
    // Esvee and Gridss
    public static final String CIPOS = "CIPOS";
    public static final String CIPOS_DESC = "Confidence interval around position for imprecise variants";

    public static final String CIRPOS = "CIRPOS";
    public static final String CIRPOS_DESC = "CIRPOS";

    public static final String SVTYPE = "SVTYPE";
    public static final String SVTYPE_DESC = "Type of structural variant";

    public static final String EVENT = "EVENT";
    public static final String EVENT_DESC = "Linking SV event ID";

    public static final String MATE_ID = "MATEID";
    public static final String MATE_ID_DESC = "Mate breakend ID";

    public static final String BEID = "BEID";
    public static final String BEID_DESC = "Identifiers of assemblies supporting the variant";

    public static final String BEIDL = "BEIDL";
    public static final String BEIDL_DESC = "Local alignment offset of corresponding BEID assembly";

    public static final String BEIDH = "BEIDL";
    public static final String BEIDH_DESC = "Remote alignment offset of corresponding BEID assembly";

    public static final String BEALN = "BEALN";
    public static final String BEALN_DESC = "Potential alignment locations of breakend sequence in the format chr:start|strand|cigar|mapq";

    public static final String HOMSEQ = "HOMSEQ";
    public static final String HOMSEQ_DESC = "Sequence of base pair identical micro-homology at event breakpoints";

    public static final String IHOMPOS = "IHOMPOS";
    public static final String IHOMPOS_DESC = "Position of inexact homology";

    public static final String QUAL = "QUAL";
    public static final String QUAL_DESC = "Quality score of breakend evidence";

    public static final String SPLIT_READS = "SR";
    public static final String SPLIT_READS_DESC = "Count of split reads supporting breakpoint";

    public static final String DISCORDANT_READS = "DP";
    public static final String DISCORDANT_READS_DESC = "Count of discordant reads supporting breakpoint";

    public static final String STRAND_BIAS = "SB";
    public static final String STRAND_BIAS_DESC = "SB";

    public static final String SV_FRAG_COUNT = "VF";
    public static final String SV_FRAG_COUNT_DESC = "Count of fragments supporting the variant breakpoint";

    public static final String ANCHOR_SUPPORT_CIGAR = "SC";
    public static final String ANCHOR_SUPPORT_CIGAR_DESC = "CIGAR of local anchor, one per assembly";

    public static final String ANCHOR_SUPPORT_CIGAR_LENGTH = "SC_LEN";
    public static final String ANCHOR_SUPPORT_CIGAR_LENGTH_DESC = "Length of local anchor CIGAR, one per assembly";


    // Esvee only
    public static final String MAPQ = "MAPQ";
    public static final String MAPQ_DESC = "The anchor's mapping quality";

    public static final String OVERHANG = "OVERHANG";
    public static final String OVERHANG_DESC = "The minimum of the left & right overhang for this assembly";

    public static final String ASSEMBLY = "ASSEMBLY";
    public static final String ASSEMBLY_DESC = "The minimum of the left & right overhang for this assembly";


    // need to check
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


    public static final String UNTEMPLATED_SEQUENCE_REPEAT_CLASS = "INSRMRC";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_TYPE = "INSRMRT";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_ORIENTATION = "INSRMRO";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_COVERAGE = "INSRMP";


    // Gridss only
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


    // Esvee only

    // also set by Gripss when run with Gridss
    public static final String LOCAL_LINKED_BY = "LOCAL_LINKED_BY";
    public static final String LOCAL_LINKED_BY_DESC = "Breakend linking information";

    public static final String REMOTE_LINKED_BY = "REMOTE_LINKED_BY";
    public static final String REMOTE_LINKED_BY_DESC = "Partner breakend linking information";

    public static final String LINKED_BY_DELIM = ",";

    // set by SvPrep depth annotation
    public static final String REF_DEPTH = "REF";
    public static final String REF_DEPTH_DESC = "Count of reads mapping across this breakend";

    public static final String REF_DEPTH_PAIR = "REFPAIR";

    public static final String REF_DEPTH_PAIR_DESC =
            "Count of reference read pairs spanning this breakend supporting the reference allele";

    public static final String ALLELE_FRACTION = "AF";
    public static final String ALLELE_FRACTION_DESC = "Allele fraction";


    // set by Gripss
    public static final String PON_FILTER_PON = "PON";
    public static final String PON_COUNT = "PON_COUNT";

    public static final String REALIGN = "REALIGN";
    public static final String REALIGN_DESC = "Variant was realigned";

    // may be merged with SV_TYPE? Gripss uses this to set StructuralVariantType
    public static final String EVENT_TYPE = "EVENTTYPE";
    public static final String EVENT_TYPE_DESC = "Structural variant type";

    public static final String ALT_PATH = "ALT_PATH";
    public static final String ALT_PATH_DESC = "Alternate path";
    public static final String RESCUE_INFO = "RESCUED";
    public static final String RESCUE_INFO_DESC = "Partner breakend rescue";
    public static final String HOTSPOT = "HOTSPOT";
    public static final String HOTSPOT_DESC = "Variant is a hotspot";

    // set by Purple
    public static final String REF_CONTEXT_FLAG = "REFG";
    public static final String REF_CONTEXT_DESC = "Reference genome surrounding break";

    public static final String RECOVERED = "RECOVERED";
    public static final String RECOVERED_DESC = "Entry has been recovered";

    public static final String RECOVERY_METHOD = "RECOVERY_METHOD";
    public static final String RECOVERY_METHOD_DESC =
            "Method used to recover, one of [UNBALANCED_SV_START, UNBALANCED_SV_END, UNSUPPORTED_BREAKEND_START, UNSUPPORTED_BREAKEND_END]";

    public static final String RECOVERY_FILTER = "RECOVERY_FILTER";
    public static final String RECOVERY_FILTER_DESC = "Filter before recovery";

    public static final String INFERRED = "INFERRED";
    public static final String INFERRED_DESC = "Breakend inferred from copy number transition";

}
