package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

class SmoothBAF {

    static double estimateBAF(@NotNull final PurpleCopyNumber region, @NotNull final PurpleCopyNumber neighbour) {
        return estimateBAF(region.averageTumorCopyNumber(), neighbour.averageActualBAF(), neighbour.averageTumorCopyNumber());
    }

    static double estimateBAF(double copyNumber, double neighbourBAF, double neighbourCopyNumber) {

        if (Doubles.isZero(copyNumber)) {
            return 0;
        }

        if (Doubles.lessOrEqual(copyNumber, 1)) {
            return 1;
        }

        double neighbourMajorAllele = neighbourBAF * neighbourCopyNumber;

        double regionMajorAllele = Math.max(0, neighbourMajorAllele - neighbourCopyNumber + copyNumber);
        double regionBAF = regionMajorAllele / copyNumber;

        return 0.5 + Math.abs(regionBAF - 0.5);
    }

}
