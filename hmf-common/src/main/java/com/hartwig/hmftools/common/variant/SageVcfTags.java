package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class SageVcfTags
{
    public static final String TIER = "TIER";
    private static final String TIER_DESCRIPTION = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";

    public static final String LOCAL_PHASE_SET = "LPS";
    private static final String PHASE_DESCRIPTION = "Local Phase Set";

    public static final String READ_CONTEXT_REPEAT_COUNT = "RC_REPC";
    private static final String READ_CONTEXT_REPEAT_COUNT_DESCRIPTION = "Repeat count at read context";

    public static final String MICROHOMOLOGY_FLAG = "MH";
    private static final String MICROHOMOLOGY_FLAG_DESCRIPTION = "Microhomology";

    public static final String TRINUCLEOTIDE_FLAG = "TNC";
    private static final String TRINUCLEOTIDE_FLAG_DESCRIPTION = "Tri-nucleotide context";

    public static final String REPEAT_COUNT_FLAG = "REP_C";
    private static final String REPEAT_COUNT_DESCRIPTION = "Repeat sequence count";

    public static final String REPEAT_SEQUENCE_FLAG = "REP_S";
    private static final String REPEAT_FLAG_DESCRIPTION = "Repeat sequence";

    public static VCFHeader addMetaData(@NotNull final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESCRIPTION));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                LOCAL_PHASE_SET, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, PHASE_DESCRIPTION));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_REPEAT_COUNT,1, VCFHeaderLineType.Integer, READ_CONTEXT_REPEAT_COUNT_DESCRIPTION));

        return header;
    }

    public static VCFHeader addRefContextHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(TRINUCLEOTIDE_FLAG, 1, VCFHeaderLineType.String, TRINUCLEOTIDE_FLAG_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_SEQUENCE_FLAG, 1, VCFHeaderLineType.String, REPEAT_FLAG_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_COUNT_FLAG, 1, VCFHeaderLineType.Integer, REPEAT_COUNT_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(MICROHOMOLOGY_FLAG, 1, VCFHeaderLineType.String, MICROHOMOLOGY_FLAG_DESCRIPTION));

        return template;
    }

}
