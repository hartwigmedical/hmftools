package com.hartwig.hmftools.common.purple.copynumber.tolerance;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

public class AlleleTolerance implements CopyNumberTolerance {

    private static final double MIN_OBSERVED_BAF_CHANGE = 0.03;
    private static final double MAX_DEVIATION_ADJUSTMENT = 0.20;
    private static final double MIN_ABSOLUTE_COPY_NUMBER_TOLERANCE = 0.3;
    private static final double MIN_RELATIVE_COPY_NUMBER_TOLERANCE = 0.1;

    @NotNull
    private final PurityAdjuster purityAdjuster;

    public AlleleTolerance(@NotNull final PurityAdjuster purityAdjuster) {
        this.purityAdjuster = purityAdjuster;
    }

    @Override
    public boolean inTolerance(@NotNull final FittedRegion first, @NotNull final FittedRegion second) {
        double purityAdjustment = purityAdjustment(purityAdjuster);

        int minBafCount = Math.min(first.bafCount(), second.bafCount());
        if (minBafCount > 0) {
            double observedBafDeviation = Math.abs(first.observedBAF() - second.observedBAF());

            double maxCopyNumber = Math.max(first.tumorCopyNumber(), second.tumorCopyNumber());
            double maxMinorAllelePloidyDeviation =
                    purityAdjustment * unadjustedMaxDeviation(MIN_ABSOLUTE_COPY_NUMBER_TOLERANCE, 0.5 * maxCopyNumber, minBafCount);
            double minorAllelePloidyDeviation = Math.abs(first.minorAllelePloidy() - second.minorAllelePloidy());
            if (Doubles.greaterThan(minorAllelePloidyDeviation, maxMinorAllelePloidyDeviation) && Doubles.greaterThan(observedBafDeviation,
                    MIN_OBSERVED_BAF_CHANGE)) {
                return false;
            }
        }

        int minWindowDepthCount = Math.min(first.depthWindowCount(), second.depthWindowCount());
        double maxDeviation = purityAdjustment * unadjustedMaxDeviation(MIN_ABSOLUTE_COPY_NUMBER_TOLERANCE, 2, minWindowDepthCount);

        boolean copyNumberInTolerance = inAbsoluteTolerance(maxDeviation, first.tumorCopyNumber(), second.tumorCopyNumber())
                || inRelativeTolerance(first.tumorCopyNumber(), second.tumorCopyNumber());

        boolean refNormalisedCopyNumberInTolerance =
                inAbsoluteTolerance(maxDeviation, first.refNormalisedCopyNumber(), second.refNormalisedCopyNumber()) || inRelativeTolerance(
                        first.refNormalisedCopyNumber(),
                        second.refNormalisedCopyNumber());

        return copyNumberInTolerance || refNormalisedCopyNumberInTolerance;
    }

    private static double purityAdjustment(@NotNull final PurityAdjuster purityAdjuster) {
        return Math.max(1, MAX_DEVIATION_ADJUSTMENT / purityAdjuster.purity());
    }

    private static double unadjustedMaxDeviation(double minTolerance, double additional, int samples) {
        return minTolerance + additional / Math.sqrt(samples);
    }

    private static boolean inAbsoluteTolerance(double tolerance, double firstCopyNumber, double secondCopyNumber) {
        final double absCopyNumberDifference = Math.abs(firstCopyNumber - secondCopyNumber);
        return Doubles.lessOrEqual(absCopyNumberDifference, tolerance);
    }

    private static boolean inRelativeTolerance(double firstCopyNumber, double secondCopyNumber) {
        return Doubles.lessOrEqual(relativeCopyNumberChange(firstCopyNumber, secondCopyNumber), MIN_RELATIVE_COPY_NUMBER_TOLERANCE);
    }

    private static double relativeCopyNumberChange(double firstCopyNumber, double secondCopyNumber) {
        final double absCopyNumberDifference = Math.abs(firstCopyNumber - secondCopyNumber);
        if (Doubles.isZero(absCopyNumberDifference)) {
            return 0;
        }

        if (Doubles.isZero(firstCopyNumber) || Doubles.isZero(secondCopyNumber)) {
            return 1;
        }

        return absCopyNumberDifference / Math.min(firstCopyNumber, secondCopyNumber);
    }

}
