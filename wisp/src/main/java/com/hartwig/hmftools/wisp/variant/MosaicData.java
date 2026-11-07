package com.hartwig.hmftools.wisp.variant;

public class MosaicData
{
    public final double TotalFragmentContribution;
    public final double TotalDualFragmentContribution;
    public final int EligibleVariantCount;

    public MosaicData(final int eligibleVariantCount, final double totalFragmentContribution, final double totalDualFragmentContribution)
    {
        EligibleVariantCount = eligibleVariantCount;
        TotalFragmentContribution = totalFragmentContribution;
        TotalDualFragmentContribution = totalDualFragmentContribution;
    }
}
