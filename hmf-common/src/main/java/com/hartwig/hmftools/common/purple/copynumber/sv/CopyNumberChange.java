package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.List;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

class CopyNumberChange {

    private final double downScale;
    private final double upScale;

    private final double downOffset;
    private final double upOffset;

    CopyNumberChange(@NotNull final List<StructuralVariantLegPloidy> structuralVariants) {
        assert (!structuralVariants.isEmpty());

        double copyNumberChangeSimple = copyNumberDifference(structuralVariants.get(0));

        int upCount = 0;
        int downCount = 0;
        double upPloidy = 0;
        double downPloidy = 0;
        double minPloidy = Double.MAX_VALUE;

        StructuralVariantLegPloidy minPloidyVariant = null;

        for (StructuralVariantLegPloidy structuralVariant : structuralVariants) {
            double ploidy = ploidy(structuralVariant);

            if (isPositive(structuralVariant)) {
                downCount++;
                downPloidy += ploidy;
            } else {
                upCount++;
                upPloidy += ploidy;
            }

            if (Doubles.lessThan(ploidy, minPloidy)) {
                minPloidy = ploidy;
                minPloidyVariant = structuralVariant;
            }
        }

        downScale = (upPloidy - copyNumberChangeSimple) / (upPloidy + downPloidy);
        upScale = 1 - downScale;

        if (upCount > 0 && downCount > 0 && minPloidyVariant != null) {
            double totalPloidy = upPloidy + downPloidy;
            double explainedPloidy = downScale * downPloidy + upScale * upPloidy;
            double offsetRequired = (totalPloidy - explainedPloidy);
            downOffset = offsetRequired / (2 * downCount);
            upOffset = downOffset * downCount / upCount;
        } else {
            downOffset = 0;
            upOffset = 0;
        }
    }

    public double copyNumberChange(@NotNull final StructuralVariantLegPloidy leg) {
        return isPositive(leg) ? downOffset + downScale * ploidy(leg) : upOffset + upScale * ploidy(leg);
    }

    private boolean isPositive(StructuralVariantLegPloidy ploidy) {
        return ploidy.orientation() == 1;
    }

    private double ploidy(StructuralVariantLegPloidy ploidy) {
        return Math.max(0, ploidy.averageImpliedPloidy());
    }

    static double copyNumberChangeSimple(@NotNull final StructuralVariantLegCopyNumber copyNumber) {

        double leftCopyNumber = copyNumber.leftCopyNumber().orElse(0D);
        double rightCopyNumber = copyNumber.rightCopyNumber().orElse(0D);

        return copyNumber.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

    private double copyNumberDifference(@NotNull final StructuralVariantLegPloidy leg) {
        return leg.rightCopyNumber().orElse(0D) - leg.leftCopyNumber().orElse(0D);
    }
}
