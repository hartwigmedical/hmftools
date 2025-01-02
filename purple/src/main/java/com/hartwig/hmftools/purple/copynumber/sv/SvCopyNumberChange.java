package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;

public class SvCopyNumberChange
{
    private final double mDownScale;
    private final double mUpScale;

    private final double mDownOffset;
    private final double mUpOffset;

    public SvCopyNumberChange(final List<StructuralVariantLegPloidy> structuralVariants)
    {
        final StructuralVariantLegPloidy template = structuralVariants.get(0);
        double copyNumberDifference = copyNumberDifference(template);

        int upCount = 0;
        int downCount = 0;
        double upPloidy = 0;
        double downPloidy = 0;
        double minPloidy = Double.MAX_VALUE;

        StructuralVariantLegPloidy minPloidyVariant = null;

        for(StructuralVariantLegPloidy structuralVariant : structuralVariants)
        {
            double ploidy = ploidy(structuralVariant);

            if(isPositive(structuralVariant))
            {
                downCount++;
                downPloidy += ploidy;
            }
            else
            {
                upCount++;
                upPloidy += ploidy;
            }

            if(Doubles.lessThan(ploidy, minPloidy))
            {
                minPloidy = ploidy;
                minPloidyVariant = structuralVariant;
            }
        }

        mDownScale = (upPloidy - copyNumberDifference) / (upPloidy + downPloidy);
        mUpScale = 1 - mDownScale;

        if(upCount > 0 && downCount > 0 && minPloidyVariant != null)
        {
            double targetDownPloidy = Math.min(template.leftCopyNumberOrZero(), downPloidy);
            double targetUpPloidy = Math.min(template.rightCopyNumberOrZero(), upPloidy);

            double targetPloidy = targetUpPloidy + targetDownPloidy;
            double explainedPloidy = mDownScale * downPloidy + mUpScale * upPloidy;
            double offsetRequired = (targetPloidy - explainedPloidy);
            mDownOffset = offsetRequired / (2 * downCount);
            mUpOffset = mDownOffset * downCount / upCount;
        }
        else
        {
            mDownOffset = 0;
            mUpOffset = 0;
        }
    }

    public double copyNumberChange(final StructuralVariantLegPloidy leg)
    {
        return isPositive(leg) ? mDownOffset + mDownScale * ploidy(leg) : mUpOffset + mUpScale * ploidy(leg);
    }

    private static boolean isPositive(final StructuralVariantLegPloidy ploidy)
    {
        return ploidy.orientation() == 1;
    }

    private static double ploidy(final StructuralVariantLegPloidy ploidy)
    {
        return Math.max(0, ploidy.averageImpliedPloidy());
    }

    static double copyNumberChangeSimple(final StructuralVariantLegCopyNumber copyNumber)
    {
        double leftCopyNumber = copyNumber.leftCopyNumberOrZero();
        double rightCopyNumber = copyNumber.rightCopyNumberOrZero();

        return copyNumber.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

    private static double copyNumberDifference(final StructuralVariantLegPloidy leg)
    {
        return leg.rightCopyNumberOrZero() - leg.leftCopyNumberOrZero();
    }
}
