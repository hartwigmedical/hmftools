package com.hartwig.hmftools.common.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;

public enum PurityAdjustment {
    ;

    public static double purityAdjustedCopyNumber(final double purity, final double normFactor, final double ratio) {
        return 2 + 2 * (ratio - normFactor) / purity / normFactor;
    }

    public static double purityAdjustedFrequency(final double purity, final double ploidy,
            final double observedFrequency, final double normalFrequency) {
        assert (greaterThan(ploidy, 0));
        assert (greaterThan(purity, 0));

        double normalPloidy = 2 * (1 - purity);
        double tumorPloidy = ploidy * purity;
        double normalAmount = 2 * (1 - purity) * normalFrequency;

        return (observedFrequency * (normalPloidy + tumorPloidy) - normalAmount) / ploidy / purity;
    }
}
