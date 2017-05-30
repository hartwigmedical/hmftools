package com.hartwig.hmftools.common.purity;

import com.hartwig.hmftools.common.numeric.Doubles;

public enum PurityAdjustment {
    ;

    public static double purityAdjustedCopynumber(double purity, double normFactor, double ratio) {
        return 2 + 2 * (ratio - normFactor) / purity / normFactor;
    }

    public static double purityAdjustedSomaticVariants(double purity, double ploidy, double observedVAF) {
        return purityAdjustedFrequency(purity, ploidy, observedVAF, 0);
    }

    public static double purityAdjustedBAF(double purity, int ploidy, double observedBAF) {
        return purityAdjustedFrequency(purity, ploidy, observedBAF, 0.5);
    }

    private static double purityAdjustedFrequency(double purity, double ploidy, double observedFrequency, double normalFrequency) {
        assert (Doubles.greaterThan(ploidy, 0));
        assert (Doubles.greaterThan(purity, 0));

        double normalPloidy = 2 * (1 - purity);
        double tumorPloidy = ploidy * purity;
        double normalAmount = 2 * (1 - purity) * normalFrequency;

        return (observedFrequency * (normalPloidy + tumorPloidy) - normalAmount) / ploidy / purity;
    }

}
