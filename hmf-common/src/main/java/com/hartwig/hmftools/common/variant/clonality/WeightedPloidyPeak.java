package com.hartwig.hmftools.common.variant.clonality;

import java.util.Collection;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public class WeightedPloidyPeak {

    private final int peakBuckets;
    private final double binWidth;


    public WeightedPloidyPeak(final int peakBuckets, final double binWidth) {
        this.peakBuckets = peakBuckets;
        this.binWidth = binWidth;
    }

    public static void doStuff(double peak, @NotNull final Collection<WeightedPloidy> ploidies) {

    }

    static double likelihoodVariantHasSpecifiedPloidy(double binWidth, double ploidy, @NotNull final WeightedPloidy weighted) {

        final BinomialDistribution binomialDistribution = new BinomialDistribution(weighted.totalReadCount(), weighted.alleleFrequency());

        double lowerBoundAlleleReadCount = Math.max(0, ploidy - binWidth / 2d) / weighted.ploidy() * weighted.alleleReadCount();
        int lowerBoundAlleleReadCountRounded = (int) Math.round(lowerBoundAlleleReadCount);
        double lowerBoundAddition = lowerBoundAlleleReadCountRounded + 0.5 - lowerBoundAlleleReadCount;

        double upperBoundAlleleReadCount = Math.max(0, ploidy + binWidth / 2d) / weighted.ploidy() * weighted.alleleReadCount();
        int upperBoundAlleleReadCountRounded = (int) Math.round(upperBoundAlleleReadCount);
        double upperBoundSubtraction = upperBoundAlleleReadCountRounded + 0.5 - upperBoundAlleleReadCount;

        return binomialDistribution.cumulativeProbability(upperBoundAlleleReadCountRounded)
                - binomialDistribution.cumulativeProbability(lowerBoundAlleleReadCountRounded)
                + lowerBoundAddition * binomialDistribution.probability(lowerBoundAlleleReadCountRounded)
                - upperBoundSubtraction * binomialDistribution.probability(upperBoundAlleleReadCountRounded);
    }

}
