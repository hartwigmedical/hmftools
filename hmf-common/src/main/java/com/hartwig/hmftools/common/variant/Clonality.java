package com.hartwig.hmftools.common.variant;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public enum Clonality {
    CLONAL,
    SUBCLONAL,
    INCONSISTENT,
    UNKNOWN;

    private static final int TRIALS = 100_000;

    public static Clonality fromSample(double copyNumber, double purity, @NotNull final AllelicDepth depth) {

        // Note: this assume normal is diploid
        double monoploidProbability = purity / (purity * copyNumber + 2 * (1 - purity));
        double monoploidSamples = depth.totalReadCount() * monoploidProbability;
        double inconsistentSamples = Math.max(copyNumber, 0) * monoploidSamples;

        final BinomialDistribution monoploidDistribution = new BinomialDistribution(TRIALS, monoploidSamples / TRIALS);
        if (depth.alleleReadCount() < monoploidDistribution.inverseCumulativeProbability(0.01)) {
            return SUBCLONAL;
        }

        final BinomialDistribution inconsistentDistribution = new BinomialDistribution(TRIALS, inconsistentSamples / TRIALS);
        if (depth.alleleReadCount() > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
            return INCONSISTENT;
        }

        return Clonality.CLONAL;
    }
}
