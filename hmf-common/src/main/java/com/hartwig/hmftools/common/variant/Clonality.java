package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public enum Clonality {
    CLONAL,
    SUBCLONAL,
    INCONSISTENT,
    UNKNOWN;

    private static final int TRIALS = 100_000;

    @NotNull
    @Deprecated
    public static Clonality fromSample(@NotNull final PurityAdjuster purityAdjuster, @NotNull final PurpleCopyNumber copyNumber,
            @NotNull final AllelicDepth depth) {
        return fromSample(purityAdjuster, copyNumber.chromosome(), copyNumber.averageTumorCopyNumber(), depth);
    }

    @NotNull
    @Deprecated
    public static Clonality fromSample(@NotNull final PurityAdjuster purityAdjuster, @NotNull final String chromosome, double copyNumber,
            @NotNull final AllelicDepth depth) {

        double purity = purityAdjuster.purity();
        int typicalCopyNumber = purityAdjuster.typicalCopyNumber(chromosome);

        double monoploidProbability = purity / (purity * copyNumber + typicalCopyNumber * (1 - purity));
        double monoploidSamples = depth.totalReadCount() * monoploidProbability;
        double inconsistentSamples = Math.max(copyNumber, 0) * monoploidSamples;

        final BinomialDistribution monoploidDistribution = new BinomialDistribution(TRIALS, Math.min(1, monoploidSamples / TRIALS));
        if (depth.alleleReadCount() < monoploidDistribution.inverseCumulativeProbability(0.01)) {
            return SUBCLONAL;
        }

        final BinomialDistribution inconsistentDistribution = new BinomialDistribution(TRIALS, Math.min(1, inconsistentSamples / TRIALS));
        if (depth.alleleReadCount() > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
            return INCONSISTENT;
        }

        return Clonality.CLONAL;
    }
}
