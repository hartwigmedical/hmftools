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
        double monoploidProbability = depth.totalReadCount() * purity / (purity * copyNumber + 2 * (1 - purity));
        double inconsistentProbability = Math.max(copyNumber, 0) * monoploidProbability;

        final BinomialDistribution monoploidDistribution = new BinomialDistribution(TRIALS, monoploidProbability / TRIALS);
        if (depth.alleleReadCount() < monoploidDistribution.inverseCumulativeProbability(0.01)) {
            return SUBCLONAL;
        }

        final BinomialDistribution inconsistentDistribution = new BinomialDistribution(TRIALS, inconsistentProbability / TRIALS);
        if (depth.alleleReadCount() > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
            return INCONSISTENT;
        }

        return Clonality.CLONAL;
    }
}
