package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

final class CopyNumberDeviation {

    @VisibleForTesting
    static final double MIN_COPY_NUMBER_TOLERANCE = 0.3;
    @VisibleForTesting
    static final double LC_MAX_COPY_NUMBER_TOLERANCE = 0.9;
    @VisibleForTesting
    static final double HC_MAX_COPY_NUMBER_TOLERANCE = 0.7;
    @VisibleForTesting
    static final int MAX_BAF_COUNT = 50;

    private static final double MAX_DEVIATION_ADJUSTMENT = 0.20;

    @NotNull
    private final PurityAdjuster purityAdjuster;

    CopyNumberDeviation(@NotNull final PurityAdjuster purityAdjuster) {
        this.purityAdjuster = purityAdjuster;
    }

    boolean inTolerance(@NotNull final FittedRegion first, @NotNull final FittedRegion second) {
        double tumorCopyNumberDeviation = Math.abs(first.tumorCopyNumber() - second.tumorCopyNumber());
        double refNormalisedCopyNumberDeviation = Math.abs(first.refNormalisedCopyNumber() - second.refNormalisedCopyNumber());
        double copyNumberDeviation = Math.min(tumorCopyNumberDeviation, refNormalisedCopyNumberDeviation);
        double rawMaxDeviation = maxCopyNumberDeviation(first, second);
        double adjustedMaxDeviation = purityAdjustedMaxCopyNumberDeviation(rawMaxDeviation);

        return Doubles.lessOrEqual(copyNumberDeviation, adjustedMaxDeviation);
    }

    static double maxCopyNumberDeviation(@NotNull final FittedRegion first, @NotNull final FittedRegion second) {
        boolean structuralBreakTransition = second.start() < first.start()
                ? !first.support().equals(SegmentSupport.NONE)
                : !second.support().equals(SegmentSupport.NONE);

        return maxCopyNumberDeviation(Math.min(first.bafCount(), second.bafCount()), structuralBreakTransition);
    }

    public double purityAdjustedMaxCopyNumberDeviation(double maxCopyNumberDeviation) {
        return maxCopyNumberDeviation * Math.max(1, MAX_DEVIATION_ADJUSTMENT / purityAdjuster.purity());
    }

    private static double maxCopyNumberDeviation(int bafCount, boolean structuralBreakTransition) {
        if (bafCount >= MAX_BAF_COUNT) {
            return MIN_COPY_NUMBER_TOLERANCE;
        }

        final double maxDeviation;
        if (structuralBreakTransition) {
            maxDeviation = HC_MAX_COPY_NUMBER_TOLERANCE;
        } else {
            maxDeviation = LC_MAX_COPY_NUMBER_TOLERANCE;
        }

        return (MIN_COPY_NUMBER_TOLERANCE - maxDeviation) / MAX_BAF_COUNT * bafCount + maxDeviation;
    }
}
