package com.hartwig.hmftools.common.sv;

public final class SvVcfTags
{
    // set by Esvee
    public static final String SVTYPE = "SVTYPE";
    public static final String SVTYPE_DESC = "Type of structural variant";

    public static final String MATE_ID = "MATEID";
    public static final String MATE_ID_DESC = "Mate breakend ID";

    public static final String CIPOS = "CIPOS";
    public static final String CIPOS_DESC = "Confidence interval around position";

    public static final String HOMSEQ = "HOMSEQ";
    public static final String HOMSEQ_DESC = "Sequence of base pair identical micro-homology at event breakpoints";

    public static final String IHOMPOS = "IHOMPOS";
    public static final String IHOMPOS_DESC = "Position of inexact homology";

    public static final String INSALN = "INSALN";
    public static final String INSALN_DESC = "Alternative alignment locations of insert sequence";

    public static final String ASMID = "ASMID";
    public static final String ASMID_DESC = "Unique id(s) of assembly(s) containing the breakend";

    public static final String ASMLEN = "ASMLEN";
    public static final String ASMLEN_DESC = "Assembly sequence length";

    public static final String ASMSEG = "ASMSEG";
    public static final String ASMSEG_DESC = "Segment indices in assembly(s) containing the breakend";

    public static final String BEAPOS = "BEAPOS";
    public static final String BEAPOS_DESC = "Breakend position(s) in assembly(s)";

    public static final String BEOR = "BEOR";
    public static final String BEOR_DESC = "Breakend orientation(s) in reference genome";

    public static final String BEAOR = "BEAOR";
    public static final String BEAOR_DESC = "Breakend orientation(s) in assembly(s)";

    public static final String SEGID = "SEGID";
    public static final String SEGID_DESC = "Unique id(s) of segment(s) containing the breakend";

    public static final String SEGALEN = "SEGALEN";
    public static final String SEGALEN_DESC = "Aligned length of segment(s) in reference genome";

    public static final String SEGMAPQ = "SEGMAPQ";
    public static final String SEGMAPQ_DESC = "MAPQ of segment containing the breakend with highest QUAL contribution";

    public static final String SEGSCO = "SEGSCO";
    public static final String SEGSCO_DESC = "Alignment score of segments containing the breakend with highest QUAL contribution";

    public static final String SEGRL = "SEGRL";
    public static final String SEGRL_DESC = "Repeat length of segment with highest QUAL contribution";

    public static final String TOTAL_FRAGS = "VF";
    public static final String TOTAL_FRAGS_DESC = "Total variant fragments supporting the breakend";

    // per sample
    public static final String SPLIT_FRAGS = "SF";
    public static final String SPLIT_FRAGS_DESC = "Count of fragments supporting the breakend with a read overlapping the breakend";

    public static final String DISC_FRAGS = "DF";
    public static final String DISC_FRAGS_DESC = "Count of discordant fragments with a read either side of the breakend";

    public static final String STRAND_BIAS = "SB";
    public static final String STRAND_BIAS_DESC = "Strand read bias";

    // set by Esvee depth annotation
    public static final String REF_DEPTH = "REF";
    public static final String REF_DEPTH_DESC = "Count of reads mapping across this breakend";

    public static final String REF_DEPTH_PAIR = "REFPAIR";

    public static final String REF_DEPTH_PAIR_DESC =
            "Count of reference read pairs spanning this breakend supporting the reference allele";

    public static final String ALLELE_FRACTION = "AF";
    public static final String ALLELE_FRACTION_DESC = "Allele frequency of the breakend";


    /*
    // also set by Gripss when run with Gridss
    public static final String LOCAL_LINKED_BY = "LOCAL_LINKED_BY";
    public static final String LOCAL_LINKED_BY_DESC = "Breakend linking information";

    public static final String REMOTE_LINKED_BY = "REMOTE_LINKED_BY";
    public static final String REMOTE_LINKED_BY_DESC = "Partner breakend linking information";

    public static final String LINKED_BY_DELIM = ",";
    */


    // set by Gripss
    public static final String PON_FILTER_PON = "PON";
    public static final String PON_COUNT = "PON_COUNT";

    public static final String HOTSPOT = "HOTSPOT";
    public static final String HOTSPOT_DESC = "Variant is a hotspot";

    public static final String UNTEMPLATED_SEQUENCE_REPEAT_CLASS = "INSRMRC";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_TYPE = "INSRMRT";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_ORIENTATION = "INSRMRO";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_COVERAGE = "INSRMP";


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
