package com.hartwig.hmftools.common.purity;

import com.hartwig.hmftools.common.numeric.Doubles;

public enum PurityAdjustment {
    ;

    public static double purityAdjustedCopyNumber(final double purity, final double normFactor, final double ratio) {
        return 2 + 2 * (ratio - normFactor) / purity / normFactor;
    }

    public static double purityAdjustedSomaticVariants(final double purity, final double ploidy,
            final double observedVAF) {
        return purityAdjustedFrequency(purity, ploidy, observedVAF, 0);
    }

    public static double purityAdjustedBAF(final double purity, final int ploidy, final double observedBAF) {
        return purityAdjustedFrequency(purity, ploidy, observedBAF, 0.5);
    }

    private static double purityAdjustedFrequency(final double purity, final double ploidy,
            final double observedFrequency, final double normalFrequency) {
        assert (Doubles.greaterThan(ploidy, 0));
        assert (Doubles.greaterThan(purity, 0));

        double normalPloidy = 2 * (1 - purity);
        double tumorPloidy = ploidy * purity;
        double normalAmount = 2 * (1 - purity) * normalFrequency;

        return (observedFrequency * (normalPloidy + tumorPloidy) - normalAmount) / ploidy / purity;
    }
}
