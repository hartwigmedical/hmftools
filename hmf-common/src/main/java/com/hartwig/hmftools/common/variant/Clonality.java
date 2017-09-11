package com.hartwig.hmftools.common.variant;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public enum Clonality {
    CLONAL,
    SUBCLONAL,
    INCONSISTENT,
    UNKNOWN;

    public static Clonality fromSample(long alleleCount, @NotNull final BinomialDistribution monoploidDistribution,
            @NotNull final BinomialDistribution inconsistentDistribution) {

        if (alleleCount < monoploidDistribution.inverseCumulativeProbability(0.01)) {
            return SUBCLONAL;
        }

        if (alleleCount > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
            return INCONSISTENT;
        }

        return Clonality.CLONAL;
    }

}
