package com.hartwig.hmftools.common.sv;

public final class SvVcfTags
{
    // set by Esvee
    public static final String ESVEE_VERSION = "EsveeVersion";

    public static final String SV_TYPE = "SVTYPE";
    public static final String SV_TYPE_DESC = "Type of structural variant";

    public static final String MATE_ID = "MATEID";
    public static final String MATE_ID_DESC = "Mate breakend ID";

    public static final String SV_ID = "SVID";
    public static final String SV_ID_DESC = "ID shared by two breakends";

    public static final String CIPOS = "CIPOS";
    public static final String CIPOS_DESC = "Confidence interval around position";

    public static final String HOMSEQ = "HOMSEQ";
    public static final String HOMSEQ_DESC = "Sequence of base pair identical micro-homology at event breakpoints";

    public static final String IHOMPOS = "IHOMPOS";
    public static final String IHOMPOS_DESC = "Position of inexact homology";

    public static final String INSALN = "INSALN";
    public static final String INSALN_DESC = "Alternative alignments of insert sequence";

    public static final String ALTALN = "ALTALN";
    public static final String ALTALN_DESC = "Alternative alignments for low map-qual breakends";

    public static final String ASM_ID = "ASMID";
    public static final String ASM_ID_DESC = "Unique id(s) of assembly(s) containing the breakend";

    public static final String ASM_LENGTH = "ASMLEN";
    public static final String ASM_LENGTH_DESC = "Assembly sequence length";

    public static final String ASM_SEG_INDEX = "ASMSEG";
    public static final String ASM_SEG_INDEX_DESC = "Segment indices in assembly(s) containing the breakend";

    public static final String BE_ASM_POS = "BEAPOS";
    public static final String BE_ASM_POS_DESC = "Breakend position(s) in assembly(s)";

    public static final String BE_ORIENT = "BEOR";
    public static final String BE_ORIENT_DESC = "Breakend orientation(s) in reference genome";

    public static final String BE_ASM_ORIENT = "BEAOR";
    public static final String BE_ASM_ORIENT_DESC = "Breakend orientation(s) in assembly(s)";

    public static final String SEG_ID = "SEGID";
    public static final String SEG_ID_DESC = "Unique id(s) of segment(s) containing the breakend";

    public static final String SEG_ALIGN_LENGTH = "SEGALEN";
    public static final String SEG_ALIGN_LENGTH_DESC = "Aligned length of segment(s) in reference genome";

    public static final String SEG_MAPQ = "SEGMAPQ";
    public static final String SEG_MAPQ_DESC = "MAPQ of segment containing the breakend with highest QUAL contribution";

    public static final String SEG_SCORE = "SEGSCO";
    public static final String SEG_SCORE_DESC = "Alignment score of segments containing the breakend with highest QUAL contribution";

    public static final String SEG_REPEAT_LENGTH = "SEGRL";
    public static final String SEG_REPEAT_LENGTH_DESC = "Repeat length of segment with highest QUAL contribution";

    // NOTE: this is used by Linx to form assembly TIs
    public static final String ASM_LINKS = "ASMLNKS";
    public static final String ASM_LINKS_DESC = "Id(s) of breakend(s) linked by assembly";

    public static final String TOTAL_FRAGS = "VF";
    public static final String TOTAL_FRAGS_DESC = "Total variant fragments supporting the breakend";

    public static final String AVG_FRAG_LENGTH = "AVGLEN";
    public static final String AVG_FRAG_LENGTH_DESC = "Average variant fragment length";

    public static final String LINE_SITE = "LINE";
    public static final String LINE_SITE_DESC = "LINE insertion site";

    public static final String UNIQUE_FRAG_POSITIONS = "UFP";
    public static final String UNIQUE_FRAG_POSITIONS_DESC = "Distinct fragment positions";

    public static final String MAX_LOCAL_REPEAT = "MLR";
    public static final String MAX_LOCAL_REPEAT_DESC = "Max local indel repeat round breakend";

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

    public static final String VCF_ITEM_DELIM = ",";

    // set by Esvee caller
    public static final String PON_FILTER_PON = "PON";
    public static final String PON_COUNT = "PON_COUNT";

    public static final String HOTSPOT = "HOTSPOT";
    public static final String HOTSPOT_DESC = "Variant is a hotspot";

    public static final String REPEAT_MASK_REPEAT_CLASS = "INSRMRC";
    public static final String REPEAT_MASK_REPEAT_CLASS_DESC = "Inserted sequence repeatmasker repeat class";
    public static final String REPEAT_MASK_REPEAT_TYPE = "INSRMRT";
    public static final String REPEAT_MASK_REPEAT_TYPE_DESC = "Inserted sequence repeatmasker repeat type";
    public static final String REPEAT_MASK_ORIENTATION = "INSRMRO";
    public static final String REPEAT_MASK_ORIENTATION_DESC = "Inserted sequence repeatmasker orientation";
    public static final String REPEAT_MASK_COVERAGE = "INSRMP";
    public static final String REPEAT_MASK_COVERAGE_DESC = "Portion of inserted sequence whose alignment overlaps the repeatmasker repeat";


    // set by Purple
    public static final String INFERRED = "INFERRED";
    public static final String INFERRED_DESC = "Breakend inferred from copy number transition";

}
