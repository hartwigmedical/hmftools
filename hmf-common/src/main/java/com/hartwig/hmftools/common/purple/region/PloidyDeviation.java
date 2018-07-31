package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.apache.commons.math3.distribution.NormalDistribution;

public class PloidyDeviation {

    private final double standardDeviation;
    private final double minStandardDevitionPerPloidyPloint;
    private final NormalDistribution dist = new NormalDistribution();

    private final double majorAlleleSubOnePenaltyMultiplier;
    private final double majorAlleleSubOneAdditionalPenalty;
    private final double baselineDeviation;

    PloidyDeviation(final double standardDeviation, final double minStandardDeviationPerPloidyPoint,
            final double majorAlleleSubOnePenaltyMultiplier, final double majorAlleleSubOneAdditionalPenalty, final double baselineDeviation) {
        this.standardDeviation = standardDeviation;
        this.minStandardDevitionPerPloidyPloint = minStandardDeviationPerPloidyPoint;
        this.majorAlleleSubOnePenaltyMultiplier = Math.abs(majorAlleleSubOnePenaltyMultiplier);
        this.majorAlleleSubOneAdditionalPenalty = Math.abs(majorAlleleSubOneAdditionalPenalty);
        this.baselineDeviation = baselineDeviation;
    }

    public double majorAlleleDeviation(final double purity, final double normFactor, final double ploidy) {
        final double majorAlleleDeviationMultiplier =
                Doubles.greaterThan(ploidy, 0) && Doubles.lessThan(ploidy, 1) ? majorAlleleSubOnePenaltyMultiplier : 1;
        final double additionalMajorAlleleDeviation =
                Math.min(majorAlleleSubOneAdditionalPenalty, Math.max(0, -majorAlleleSubOneAdditionalPenalty * (ploidy - 1)));
        return majorAlleleDeviationMultiplier * alleleDeviation(purity, normFactor, ploidy) + additionalMajorAlleleDeviation;
    }

    public double minorAlleleDeviation(final double purity, final double normFactor, final double ploidy) {
        return alleleDeviation(purity, normFactor, ploidy);
    }

    private double alleleDeviation(final double purity, final double normFactor, final double ploidy) {
        if (Doubles.lessThan(ploidy, -0.5)) {
            return 1 + baselineDeviation;
        }

        double ploidyDistanceFromInteger = Doubles.absDistanceFromInteger(ploidy);
        double standardDeviationsPerPloidy = Math.max(minStandardDevitionPerPloidyPloint, purity * normFactor / 2 / standardDeviation);
        return 2 * dist.cumulativeProbability(ploidyDistanceFromInteger * standardDeviationsPerPloidy) - 1 + baselineDeviation;
    }

}
