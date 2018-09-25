package com.hartwig.hmftools.common.purple.copynumber.tolerance;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

public class AlleleDeviation implements CopyNumberTolerance {

    @VisibleForTesting
    static final double MIN_COPY_NUMBER_TOLERANCE = 0.3;
    @VisibleForTesting
    static final double MAX_COPY_NUMBER_TOLERANCE = 0.9;
    @VisibleForTesting
    static final double MAX_DEVIATION_ADJUSTMENT = 0.20;
    @VisibleForTesting
    static final int MAX_BAF_COUNT = 50;

    @NotNull
    private final PurityAdjuster purityAdjuster;

    AlleleDeviation(@NotNull final PurityAdjuster purityAdjuster) {
        this.purityAdjuster = purityAdjuster;
    }

    @Override
    public boolean inTolerance(@NotNull final FittedRegion first, @NotNull final FittedRegion second) {
        int minBafCount = Math.min(first.bafCount(), second.bafCount());
        double maxDeviation = purityAdjustedMaxDeviation(purityAdjuster, minBafCount);

        if (minBafCount > 0) {
            double minorAllelePloidyDeviation = Math.abs(first.minorAllelePloidy() - second.minorAllelePloidy());
            if (Doubles.greaterThan(minorAllelePloidyDeviation, maxDeviation)) {
                return false;
            }
        }

        double tumorCopyNumberDeviation = Math.abs(first.tumorCopyNumber() - second.tumorCopyNumber());
        double refNormalisedCopyNumberDeviation = Math.abs(first.refNormalisedCopyNumber() - second.refNormalisedCopyNumber());
        double copyNumberDeviation = Math.min(tumorCopyNumberDeviation, refNormalisedCopyNumberDeviation);

        return Doubles.lessOrEqual(copyNumberDeviation, maxDeviation);
    }

    private static double purityAdjustedMaxDeviation(@NotNull final PurityAdjuster purityAdjuster, int bafCount) {
        double rawMaxDeviation = unadjustedMaxDeviation(bafCount);
        return rawMaxDeviation * Math.max(1, MAX_DEVIATION_ADJUSTMENT / purityAdjuster.purity());
    }

    private static double unadjustedMaxDeviation(int bafCount) {
        if (bafCount >= MAX_BAF_COUNT) {
            return MIN_COPY_NUMBER_TOLERANCE;
        }

        return (MIN_COPY_NUMBER_TOLERANCE - MAX_COPY_NUMBER_TOLERANCE) / MAX_BAF_COUNT * bafCount + MAX_COPY_NUMBER_TOLERANCE;
    }

}
