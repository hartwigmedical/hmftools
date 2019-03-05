package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class CopyNumberChange {

    private final double totalDownPloidy;
    private final double totalUpPloidy;

    private final double positiveCopyNumberProportion;
    private final double negativeCopyNumberProportion;

    CopyNumberChange(@NotNull final List<StructuralVariantLegPloidy> structuralVariants) {

//        if (structuralVariants.stream().anyMatch(x -> x.position() == 103821416L)) {
//            System.out.println("sdf");
//        }

        totalDownPloidy = structuralVariants.stream().filter(this::isPositive).mapToDouble(this::ploidy).sum();
        totalUpPloidy = structuralVariants.stream().filter(this::isNegative).mapToDouble(this::ploidy).sum();
        final double totalPloidy = Math.abs(totalDownPloidy - totalUpPloidy);
        final double maxPloidy = Math.max(totalDownPloidy, totalUpPloidy);

        positiveCopyNumberProportion = totalDownPloidy / totalPloidy;
        negativeCopyNumberProportion = totalUpPloidy / totalPloidy;
    }

    public double copyNumberChange(@NotNull final StructuralVariantLegPloidy leg) {
//        double copyNumberChangeSimple = Math.abs(copyNumberChangeSimple(leg));
        double copyNumberChangeSimple = leg.rightCopyNumber().orElse(0D) - leg.leftCopyNumber().orElse(0D);

        double downScale = (totalUpPloidy - copyNumberChangeSimple) / (totalUpPloidy + totalDownPloidy);
        double upScale = 1 - downScale;

        return isPositive(leg)
                ? 0.25+ downScale * ploidy(leg)
                : 0.5 + upScale * ploidy(leg);

//
//        return isPositive(leg)
//                ? Math.abs(copyNumberChangeSimple) * positiveCopyNumberProportion * ploidy(leg) / totalPositivePloidy
//                : Math.abs(copyNumberChangeSimple) * negativeCopyNumberProportion * ploidy(leg) / totalNegativePloidy;
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
