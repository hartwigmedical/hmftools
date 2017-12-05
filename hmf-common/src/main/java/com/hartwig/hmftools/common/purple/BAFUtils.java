package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.numeric.Doubles;

public class BAFUtils {

    public static final double NORMAL_BAF = 0.535;

    public static int minAlleleCount(final int ploidy) {
        return (int) Math.max(0, Math.round(ploidy / 2d));
    }

    public static double modelBAF(final double tumorPurity, final int ploidy, final int betaAlleleCount) {
        assert (betaAlleleCount >= ploidy / 2d);
        double normalPurity = 1 - tumorPurity;

        double betaObservations = 1 * normalPurity + betaAlleleCount * tumorPurity;
        double totalObservations = 2 * normalPurity + ploidy * tumorPurity;

        if (Doubles.isZero(totalObservations)) {
            return NORMAL_BAF;
        }

        return Math.max(NORMAL_BAF, betaObservations / totalObservations);
    }
}
