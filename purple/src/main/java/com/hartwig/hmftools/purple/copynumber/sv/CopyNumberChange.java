package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class CopyNumberChange
{
    private final double mDownScale;
    private final double mUpScale;

    private final double mDownOffset;
    private final double mUpOffset;

    public CopyNumberChange(@NotNull final List<StructuralVariantLegPloidy> structuralVariants)
    {
        assert (!structuralVariants.isEmpty());

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
            double targetDownPloidy = Math.min(template.leftCopyNumber().orElse(0d), downPloidy);
            double targetUpPloidy = Math.min(template.rightCopyNumber().orElse(0d), upPloidy);

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

    public double copyNumberChange(@NotNull final StructuralVariantLegPloidy leg)
    {
        return isPositive(leg) ? mDownOffset + mDownScale * ploidy(leg) : mUpOffset + mUpScale * ploidy(leg);
    }

    private static boolean isPositive(StructuralVariantLegPloidy ploidy)
    {
        return ploidy.orientation() == 1;
    }

    private static double ploidy(StructuralVariantLegPloidy ploidy)
    {
        return Math.max(0, ploidy.averageImpliedPloidy());
    }

    static double copyNumberChangeSimple(@NotNull final StructuralVariantLegCopyNumber copyNumber)
    {
        double leftCopyNumber = copyNumber.leftCopyNumber().orElse(0D);
        double rightCopyNumber = copyNumber.rightCopyNumber().orElse(0D);

        return copyNumber.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

    private static double copyNumberDifference(@NotNull final StructuralVariantLegPloidy leg)
    {
        return leg.rightCopyNumber().orElse(0D) - leg.leftCopyNumber().orElse(0D);
    }
}
