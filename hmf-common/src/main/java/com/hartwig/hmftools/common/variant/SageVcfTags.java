package com.hartwig.hmftools.common.variant;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class SageVcfTags
{
    public static final String TIER = "TIER";
    public static final String TIER_DESC = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";

    public static final String LOCAL_PHASE_SET = "LPS";
    public static final String LOCAL_PHASE_SET_DESC = "Local Phase Set";

    public static final String READ_CONTEXT_REPEAT_COUNT = "RC_REPC";
    public static final String READ_CONTEXT_REPEAT_COUNT_DESC = "Repeat count at read context";

    public static final String MICROHOMOLOGY_FLAG = "MH";
    public static final String MICROHOMOLOGY_FLAG_DESCRIPTION = "Microhomology";

    public static final String TRINUCLEOTIDE_FLAG = "TNC";
    public static final String TRINUCLEOTIDE_FLAG_DESCRIPTION = "Tri-nucleotide context";

    public static final String REPEAT_COUNT_FLAG = "REP_C";
    public static final String REPEAT_COUNT_DESCRIPTION = "Repeat sequence count";

    public static final String REPEAT_SEQUENCE_FLAG = "REP_S";
    public static final String REPEAT_FLAG_DESCRIPTION = "Repeat sequence";

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

    public static VCFHeader addRefContextHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(TRINUCLEOTIDE_FLAG, 1, VCFHeaderLineType.String, TRINUCLEOTIDE_FLAG_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_SEQUENCE_FLAG, 1, VCFHeaderLineType.String, REPEAT_FLAG_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_COUNT_FLAG, 1, VCFHeaderLineType.Integer, REPEAT_COUNT_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(MICROHOMOLOGY_FLAG, 1, VCFHeaderLineType.String, MICROHOMOLOGY_FLAG_DESCRIPTION));

        return template;
    }
}
