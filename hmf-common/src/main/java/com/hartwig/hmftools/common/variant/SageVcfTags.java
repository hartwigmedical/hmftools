package com.hartwig.hmftools.common.variant;

import static java.lang.String.format;

import com.hartwig.hmftools.common.bam.ConsensusType;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public final class SageVcfTags
{
    public static final String TIER = "TIER";
    public static final String TIER_DESC = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";

    public static final String LOCAL_PHASE_SET = "LPS";
    public static final String LOCAL_PHASE_SET_DESC = "Local Phase Set";

    public static final String LPS_APPEND_INFO = "LPSA";
    public static final String LPS_APPEND_INFO_DESC = "Local Phase Set Append Info";

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

    public static final String MAP_QUAL_FACTOR = "MQF";
    public static final String MAP_QUAL_FACTOR_DESC = "Map qual heuristic as used in min tumor quality filter";

    public static final String AVG_RAW_BASE_QUAL = "RABQ";
    public static final String AVG_RAW_BASE_QUAL_DESC = "Average calculated raw base quality in alt reads";

    public static final String AVG_BASE_QUAL = "ABQ";
    public static final String AVG_BASE_QUAL_DESC = "Average calculated base quality (all,alt)";

    public static final String NEARBY_INDEL_FLAG = "NEARBY_INDEL";
    public static final String NEARBY_INDEL_FLAG_DESC = "Variant has an INDEL overlapping its core";

    public static final String MIN_COORDS_COUNT = "MUC";
    public static final String MIN_COORDS_COUNT_DESC = "Min unique fragment coordinates in alt reads";

    // deprecated after v4.1, replaced by percentage
    public static final String AVG_READ_EDGE_DISTANCE = "AED";
    public static final String AVG_READ_EDGE_DISTANCE_DESC = "Average read edge distance [alt,total]";

    public static final String AVG_EDGE_DISTANCE_PERC = "AED";
    public static final String AVG_EDGE_DISTANCE_PERC_DESC = "Average edge distance as percent of read length [alt,total]";

    public static final int CONSENSUS_TYPE_COUNT = ConsensusType.values().length;
    public static final int CONSENSUS_TAG_TYPE_COUNT = CONSENSUS_TYPE_COUNT * 2;

    public static final String TINC_LEVEL = "tincLevel";

    public static final String TINC_RECOVERED_FLAG = "TINC_RECOVERED";
    public static final String TINC_RECOVERED_DESC = "Variant recovered from germline filters by TINC detection";

    public static final String LIST_SEPARATOR = ",";

    public static void writeTincLevel(final VCFHeader vcfHeader, final double tincLevel)
    {
        if(tincLevel > 0)
            vcfHeader.addMetaDataLine(new VCFHeaderLine(TINC_LEVEL, format("%.3f", tincLevel)));
    }

    public static double parseTincLevel(final VCFHeader vcfHeader)
    {
        VCFHeaderLine tincHeader = vcfHeader.getMetaDataLine(TINC_LEVEL);

        if(tincHeader == null)
            return 0;

        try
        {
            return Double.parseDouble(tincHeader.getValue());
        }
        catch(Exception e)
        {
            return 0;
        }
    }
}
