package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.apache.commons.math3.distribution.NormalDistribution;

public class PloidyDeviation {

    private final double standardDeviation;
    private final NormalDistribution dist = new NormalDistribution();

    PloidyDeviation(final double standardDeviation) {
        this.standardDeviation = standardDeviation;
    }

    public double majorAlleleDeivation(final double purity, final double normFactor, final double ploidy) {
        return Doubles.lessThan(ploidy, 0.5) ? 1 : minorAlleleDeviation(purity, normFactor, ploidy);
    }

    public double minorAlleleDeviation(final double purity, final double normFactor, final double ploidy) {
        double probability = Math.abs(Math.round(ploidy) - ploidy) * purity * normFactor / 2 / standardDeviation;
        return 2 * dist.cumulativeProbability(probability) - 1;
    }

}
