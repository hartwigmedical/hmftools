package com.hartwig.hmftools.common.purple;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.baf.ExpectedBAF;

public class BAFUtils {

    @Deprecated
    public static final double NORMAL_BAF = 0.535;

    private static final double AMBIGUOUS_OFFSET = 0.7;
    private final double expectedBaf;

    public BAFUtils(final double expectedBaf) {
        this.expectedBaf = expectedBaf;
    }

    public BAFUtils(final int averageDepth) {
        this.expectedBaf = ExpectedBAF.expectedBAF(averageDepth);
    }

    public double expectedBAF() {
        return expectedBaf;
    }

    public double ambiguousBAF() {
        return expectedBaf + AMBIGUOUS_OFFSET;
    }

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

    @VisibleForTesting
    public static double[] modelBAFToMinimizeDeviation(final double purity, final int ploidy, final double actualBAF) {
        double finalBAF = 0;
        double finalDeviation = 0;
        int finalMajorAllele = 0;

        int minBetaAllele = BAFUtils.minAlleleCount(ploidy);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            double modelBAF = BAFUtils.modelBAF(purity, ploidy, betaAllele);
            double modelDeviation = bafDeviation(modelBAF, actualBAF);

            if (betaAllele == minBetaAllele || modelDeviation < finalDeviation) {
                finalBAF = modelBAF;
                finalDeviation = modelDeviation;
                finalMajorAllele = betaAllele;
            }
        }

        return new double[] { finalBAF, finalDeviation, finalMajorAllele };
    }


    @VisibleForTesting
    public static double bafDeviation(final double modelBAF, final double actualBAF) {
        return Math.abs(modelBAF - actualBAF);
    }

    private boolean isAmbiguous(final double frequency) {
        return Doubles.lessThan(frequency, ambiguousBAF());
    }

}
