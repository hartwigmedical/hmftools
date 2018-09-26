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
    static final double MAX_DEVIATION_ADJUSTMENT = 0.20;

    @NotNull
    private final PurityAdjuster purityAdjuster;

    AlleleDeviation(@NotNull final PurityAdjuster purityAdjuster) {
        this.purityAdjuster = purityAdjuster;
    }

    @Override
    public boolean inTolerance(@NotNull final FittedRegion first, @NotNull final FittedRegion second) {
        double purityAdjustment = purityAdjustment(purityAdjuster);

        int minBafCount = Math.min(first.bafCount(), second.bafCount());
        if (minBafCount > 0) {
            double maxCopyNumber = Math.max(first.tumorCopyNumber(), second.tumorCopyNumber());
            double maxBafDeviation = purityAdjustment * unadjustedMaxDeviation(MIN_COPY_NUMBER_TOLERANCE, 0.5 * maxCopyNumber, minBafCount);
            double minorAllelePloidyDeviation = Math.abs(first.minorAllelePloidy() - second.minorAllelePloidy());
            if (Doubles.greaterThan(minorAllelePloidyDeviation, maxBafDeviation)) {
                return false;
            }
        }

        int minWindowDepthCount = Math.min(first.depthWindowCount(), second.depthWindowCount());
        double maxDeviation = purityAdjustment * unadjustedMaxDeviation(MIN_COPY_NUMBER_TOLERANCE, 2, minWindowDepthCount);
        double tumorCopyNumberDeviation = Math.abs(first.tumorCopyNumber() - second.tumorCopyNumber());
        double refNormalisedCopyNumberDeviation = Math.abs(first.refNormalisedCopyNumber() - second.refNormalisedCopyNumber());
        double copyNumberDeviation = Math.min(tumorCopyNumberDeviation, refNormalisedCopyNumberDeviation);

        return Doubles.lessOrEqual(copyNumberDeviation, maxDeviation);
    }

    private static double purityAdjustment(@NotNull final PurityAdjuster purityAdjuster) {
        return Math.max(1, MAX_DEVIATION_ADJUSTMENT / purityAdjuster.purity());
    }

    private static double unadjustedMaxDeviation(double minTolerance, double additional, int samples) {
        return minTolerance + additional / Math.sqrt(samples);
    }

}
