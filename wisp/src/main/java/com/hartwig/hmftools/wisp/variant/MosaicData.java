package com.hartwig.hmftools.wisp.variant;

public class MosaicData
{
    public final double TotalFragmentContribution;
    public final double TotalDualFragmentContribution;
    public final int EligibleVariantCount;

    public static final MosaicData INVALID_DATA = new MosaicData(0, 0, 0);

    public MosaicData(final int eligibleVariantCount, final double totalFragmentContribution, final double totalDualFragmentContribution)
    {
        EligibleVariantCount = eligibleVariantCount;
        TotalFragmentContribution = totalFragmentContribution;
        TotalDualFragmentContribution = totalDualFragmentContribution;
    }

    public double fragmentsPerEligibleVariant()
    {
        return EligibleVariantCount > 0 ? TotalFragmentContribution / EligibleVariantCount : 0;
    }

    public double dualFragmentsPerEligibleVariant()
    {
        return EligibleVariantCount > 0 ? TotalDualFragmentContribution / EligibleVariantCount : 0;
    }
}
