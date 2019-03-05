package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class CopyNumberChange {

    private final double totalPositivePloidy;
    private final double totalNegativePloidy;

    private final double positiveCopyNumberProportion;
    private final double negativeCopyNumberProportion;

    CopyNumberChange(@NotNull final List<StructuralVariantLegPloidy> structuralVariants) {

        totalPositivePloidy = structuralVariants.stream().filter(this::isPositive).mapToDouble(this::ploidy).sum();
        totalNegativePloidy = structuralVariants.stream().filter(this::isNegative).mapToDouble(this::ploidy).sum();
        final double totalPloidy = totalPositivePloidy + totalNegativePloidy;

        positiveCopyNumberProportion = totalPositivePloidy / totalPloidy;
        negativeCopyNumberProportion = totalNegativePloidy / totalPloidy;
    }

    public double copyNumberChange(@NotNull final StructuralVariantLegPloidy leg) {
        double copyNumberChangeSimple = copyNumberChangeSimple(leg);

        return isPositive(leg)
                ? copyNumberChangeSimple * positiveCopyNumberProportion * ploidy(leg) / totalPositivePloidy
                : copyNumberChangeSimple * negativeCopyNumberProportion * ploidy(leg) / totalNegativePloidy;
    }

    private boolean isPositive(StructuralVariantLegPloidy ploidy) {
        return ploidy.orientation() == 1;
    }

    private boolean isNegative(StructuralVariantLegPloidy ploidy) {
        return !isPositive(ploidy);
    }

    private double ploidy(StructuralVariantLegPloidy ploidy) {
        return Math.max(0, ploidy.averageImpliedPloidy());
    }

    static double copyNumberChangeSimple(@NotNull final StructuralVariantLegCopyNumber copyNumber) {

        double leftCopyNumber = copyNumber.leftCopyNumber().orElse(0D);
        double rightCopyNumber = copyNumber.rightCopyNumber().orElse(0D);

        return copyNumber.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }
}
