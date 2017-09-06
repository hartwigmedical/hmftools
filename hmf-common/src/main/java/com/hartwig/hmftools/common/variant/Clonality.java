package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public enum Clonality {
    CLONAL,
    SUBCLONAL,
    INCONSISTENT,
    UNKNOWN;

    private static final int MULTIPLIER = 1_000;

    public static Clonality fromSample(double purity, double normFactor, double observedNormalRatio, double observedTumorRatio,
            @NotNull AllelicDepth depth) {

        try {
            if (Doubles.isZero(observedTumorRatio)) {
                return UNKNOWN;
            }

            double monoploidProbability = purity * normFactor / 2 / observedTumorRatio;
            final BinomialDistribution monoploidDistribution =
                    new BinomialDistribution(depth.totalReadCount() * MULTIPLIER, monoploidProbability / MULTIPLIER);
            if (depth.alleleReadCount() < monoploidDistribution.inverseCumulativeProbability(0.01)) {
                return SUBCLONAL;
            }

            double inconsistentProbability = (purity - 1) * observedNormalRatio * normFactor / observedTumorRatio + 1;
            final BinomialDistribution inconsistentDistribution =
                    new BinomialDistribution(depth.totalReadCount() * MULTIPLIER, inconsistentProbability / MULTIPLIER);
            if (depth.alleleReadCount() > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
                return INCONSISTENT;
            }

            return CLONAL;
        } catch (Exception e) {
            return UNKNOWN;
        }
    }
}
