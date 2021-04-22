package com.hartwig.hmftools.purple.region;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.NormalDistribution;

class PloidyDeviation {

    private final double standardDeviation;
    private final double minStandardDeviationPerPloidyPoint;
    private final NormalDistribution dist = new NormalDistribution();

    private final double majorAlleleSubOnePenaltyMultiplier;
    private final double majorAlleleSubMinAdditionalPenalty;
    private final double minDeviation;

    PloidyDeviation(final double standardDeviation, final double minStandardDeviationPerPloidyPoint,
            final double majorAlleleSubOnePenaltyMultiplier, final double majorAlleleSubOneAdditionalPenalty, final double minDeviation) {
        this.standardDeviation = standardDeviation;
        this.minStandardDeviationPerPloidyPoint = minStandardDeviationPerPloidyPoint;
        this.majorAlleleSubOnePenaltyMultiplier = Math.abs(majorAlleleSubOnePenaltyMultiplier);
        this.majorAlleleSubMinAdditionalPenalty = Math.abs(majorAlleleSubOneAdditionalPenalty);
        this.minDeviation = minDeviation;
    }

    double majorAlleleDeviation(final double purity, final double normFactor, final double ploidy) {
        final double majorAlleleDeviationMultiplier = Doubles.greaterThan(ploidy, 0) && Doubles.lessThan(ploidy, 1)
                ? Math.max(1, majorAlleleSubOnePenaltyMultiplier * (1 - ploidy))
                : 1;
        final double deviation =
                majorAlleleDeviationMultiplier * alleleDeviation(purity, normFactor, ploidy) + subMinAdditionalPenalty(1, ploidy);
        return Math.max(deviation, minDeviation);
    }

    double minorAlleleDeviation(final double purity, final double normFactor, final double ploidy) {
        final double deviation = alleleDeviation(purity, normFactor, ploidy) + subMinAdditionalPenalty(0, ploidy);
        return Math.max(deviation, minDeviation);
    }

    private double alleleDeviation(final double purity, final double normFactor, final double ploidy) {
        final double ploidyDistanceFromInteger = Doubles.lessThan(ploidy, -0.5)
                ? 0.5
                : Doubles.absDistanceFromInteger(ploidy);

        double standardDeviationsPerPloidy = Math.max(minStandardDeviationPerPloidyPoint, purity * normFactor / 2 / standardDeviation);
        return 2 * dist.cumulativeProbability(ploidyDistanceFromInteger * standardDeviationsPerPloidy) - 1 + Math.max(-0.5 - ploidy, 0);
    }

    private double subMinAdditionalPenalty(final double minPloidy, final double ploidy) {
        return Math.min(majorAlleleSubMinAdditionalPenalty, Math.max(0, -majorAlleleSubMinAdditionalPenalty * (ploidy - minPloidy)));
    }
}
