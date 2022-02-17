package com.hartwig.hmftools.common.sage;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class SageMetaData
{
    public static final String TIER = "TIER";
    public static final String LOCAL_PHASE_SET = "LPS";

    private static final String TIER_DESCRIPTION = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";
    private static final String PHASE_DESCRIPTION = "Local Phase Set";
    private static final String REALIGN_DESCRIPTION = "Local Realignment Set";

    private SageMetaData() { }

    @NotNull
    public static VCFHeader addSageMetaData(@NotNull final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(LOCAL_PHASE_SET, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, PHASE_DESCRIPTION));
        return header;
    }
}
