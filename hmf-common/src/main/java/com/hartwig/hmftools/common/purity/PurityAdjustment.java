package com.hartwig.hmftools.common.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;

public enum PurityAdjustment {
    ;

    public static double purityAdjustedCopyNumber(final Gender gender, final String chromosome, final double purity, final double normFactor,
            final double ratio) {
        final double expectedCopyNumber = gender.equals(Gender.MALE) && (chromosome.equals("X") || chromosome.equals("Y")) ? 1 : 2;
        return purityAdjustedCopyNumber(purity, normFactor, ratio, expectedCopyNumber);
    }

    public static double purityAdjustedCopyNumber(final double purity, final double normFactor, final double ratio, final double expectedCopyNumber) {
        return Doubles.isZero(ratio) ? 0 : expectedCopyNumber + (2 * ratio - expectedCopyNumber * normFactor) / purity / normFactor;
    }

    public static double purityAdjustedVAF(final double purity, final double copyNumber, final double observedFrequency) {
        return purityAdjustedFrequency(purity, copyNumber, observedFrequency, 0);
    }

    public static double purityAdjustedBAF(final double purity, final double copyNumber, final double observedFrequency) {
        return purityAdjustedFrequency(purity, copyNumber, observedFrequency, 0.5);
    }

    private static double purityAdjustedFrequency(final double purity, final double copyNumber, final double observedFrequency, final double normalFrequency) {
        assert (greaterThan(copyNumber, 0));
        assert (greaterThan(purity, 0));

        double normalPloidy = 2 * (1 - purity);
        double tumorPloidy = copyNumber * purity;
        double normalAmount = 2 * (1 - purity) * normalFrequency;

        return (observedFrequency * (normalPloidy + tumorPloidy) - normalAmount) / copyNumber / purity;
    }
}
