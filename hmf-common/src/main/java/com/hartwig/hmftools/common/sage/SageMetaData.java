package com.hartwig.hmftools.common.sage;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class SageMetaData
{
    public static final String TIER = "TIER";
    public static final String LOCAL_PHASE_SET = "LPS";
    public static final String LOCAL_REALIGN_SET = "LRS";

    private static final String TIER_DESCRIPTION = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";
    private static final String PHASE_DESCRIPTION = "Local Phase Set";
    private static final String REALIGN_DESCRIPTION = "Local Realignment Set";

    // no longer required by Pave nor set by Sage, to be deprecated once SnpEff handling is removed
    public static final String PHASED_INFRAME_INDEL = "PII";
    private static final String PHASED_INFRAME_INDEL_DESCRIPTION = "Phased inframe indel identifier";

    private SageMetaData()
    {
    }

    @NotNull
    public static VCFHeader addSageMetaData(@NotNull final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(LOCAL_PHASE_SET, 1, VCFHeaderLineType.Integer, PHASE_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(LOCAL_REALIGN_SET, 1, VCFHeaderLineType.Integer, REALIGN_DESCRIPTION));
        return header;
    }
}
