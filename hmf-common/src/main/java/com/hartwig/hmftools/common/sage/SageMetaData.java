package com.hartwig.hmftools.common.sage;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SageMetaData {

    public final static String TIER = "TIER";
    public final static String PHASED_INFRAME_INDEL = "PII";


    private final static String TIER_DESCRIPTION = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";
    public final static String PHASED_INFRAME_INDEL_DESCRIPTION = "Phased inframe indel identifier";

    @NotNull
    public static VCFHeader addSageMetaData(@NotNull final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(PHASED_INFRAME_INDEL, 1, VCFHeaderLineType.Integer, PHASED_INFRAME_INDEL_DESCRIPTION));
        return header;
    }

}
