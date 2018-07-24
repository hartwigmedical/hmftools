package com.hartwig.hmftools.common.purple.region;

import org.apache.commons.math3.distribution.NormalDistribution;

public class PloidyDeviation {

    private final double purity;
    private final double normFactor;
    private final double standardDeviation;

    public PloidyDeviation(final double purity, final double normFactor, final double standardDeviation) {
        this.purity = purity;
        this.normFactor = normFactor;
        this.standardDeviation = standardDeviation;
    }

    public double deviation(double ploidy) {
        double probability = Math.abs(Math.round(ploidy) - ploidy) * purity * normFactor / 2;
        NormalDistribution dist = new NormalDistribution(0, standardDeviation);
        return 2 * dist.cumulativeProbability(probability) - 1;
    }

}
