package com.hartwig.hmftools.purple.copynumber.sv;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

public class StructuralVariantLegPloidy extends StructuralVariantLegCopyNumber
{
    private final double mObservedVaf;
    private final double mAdjustedVaf;
    private final double mUnweightedImpliedPloidy;

    private double mAverageImpliedPloidy;
    private double mWeight;

    public StructuralVariantLegPloidy(
            final StructuralVariantLegCopyNumber leg, final double observedVaf, final double adjustedVaf,
            final double averageImpliedPloidy, final double unweightedImpliedPloidy, final double weight)
    {
        super(leg, leg.leftCopyNumber(), leg.rightCopyNumber());
        mObservedVaf = observedVaf;
        mAdjustedVaf = adjustedVaf;
        mAverageImpliedPloidy = averageImpliedPloidy;
        mUnweightedImpliedPloidy = unweightedImpliedPloidy;
        mWeight = weight;
    }

    public double observedVaf() { return mObservedVaf; }
    public double adjustedVaf() { return mAdjustedVaf; }
    public double unweightedImpliedPloidy() { return mUnweightedImpliedPloidy; }

    public double averageImpliedPloidy() { return mAverageImpliedPloidy; }
    public void setAverageImpliedPloidy(double ploidy) { mAverageImpliedPloidy = ploidy; }

    public double weight() { return mWeight; }
    public void setWeight(double weight) { mWeight = weight; }

    public double impliedRightCopyNumberWeight()
    {
        return canInferRight() ? weight() : 0;
    }

    public double impliedRightCopyNumber()
    {
        if(isDecreasingFromZero(1, leftCopyNumber()))
            return 0;

        return leftCopyNumber() != null ? leftCopyNumber() - orientation() * averageImpliedPloidy() : 0;
    }

    public double impliedLeftCopyNumber()
    {
        if(isDecreasingFromZero(-1, rightCopyNumber()))
            return 0;

        return rightCopyNumber() != null ? rightCopyNumber() + orientation() * averageImpliedPloidy() : 0;
    }

    public double impliedLeftCopyNumberWeight()
    {
        return canInferLeft() ? weight() : 0;
    }

    private boolean canInferRight()
    {
        return isDecreasingFromZero(1, leftCopyNumber()) || leftCopyNumber() != null;
    }

    private boolean canInferLeft()
    {
        return isDecreasingFromZero(-1, rightCopyNumber()) || rightCopyNumber() != null;
    }

    private boolean isDecreasingFromZero(int orientation, final Double copyNumber)
    {
        return orientation() == orientation && copyNumber != null && Doubles.isZero(copyNumber);
    }
}
