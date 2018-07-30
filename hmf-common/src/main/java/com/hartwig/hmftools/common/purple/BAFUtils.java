package com.hartwig.hmftools.common.purple;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.baf.ExpectedBAF;

import org.jetbrains.annotations.NotNull;

public class BAFUtils {

    private static final double AMBIGUOUS_OFFSET = 0.007;
    private static final double CLONAL_DISTANCE = 0.25;
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

    public double modelBAF(final double tumorPurity, final int ploidy, final int betaAlleleCount) {
        assert (betaAlleleCount >= ploidy / 2d);
        double normalPurity = 1 - tumorPurity;

        double betaObservations = 1 * normalPurity + betaAlleleCount * tumorPurity;
        double totalObservations = 2 * normalPurity + ploidy * tumorPurity;

        if (Doubles.isZero(totalObservations)) {
            return expectedBaf;
        }

        return Math.max(expectedBaf, betaObservations / totalObservations);
    }

    @VisibleForTesting
    public double[] modelBAFToMinimizeDeviation(final double purity, final int ploidy, final double actualBAF) {
        double finalBAF = 0;
        double finalDeviation = 0;
        int finalMajorAllele = 0;

        int minBetaAllele = BAFUtils.minAlleleCount(ploidy);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            double modelBAF = modelBAF(purity, ploidy, betaAllele);
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

    public double purityAdjustedBAF(@NotNull PurityAdjuster adjuster, @NotNull final String chromosomeName, final double copyNumber,
            final double observedFrequency) {
        int normalCopyNumber = adjuster.typicalCopyNumber(chromosomeName);
        if (normalCopyNumber == 1 || (Doubles.positive(observedFrequency) && Doubles.lessOrEqual(copyNumber, 1))) {
            return 1;
        }

        if (Doubles.isZero(observedFrequency)) {
            return 0;
        }

        double typicalFrequency = 0.5;
        double rawAdjustedBaf = adjuster.purityAdjustedFrequency(copyNumber, observedFrequency, normalCopyNumber, typicalFrequency);

        int ploidy = (int) Math.round(copyNumber);
        if (isAmbiguous(observedFrequency) && isClonal(copyNumber) && ploidy > 0) {
            int minBetaAllele = BAFUtils.minAlleleCount(ploidy);
            double modelBAF = modelBAF(adjuster.purity(), ploidy, minBetaAllele);
            if (isAmbiguous(modelBAF)) {
                return (double) minBetaAllele / ploidy;
            }
        }

        return rawAdjustedBaf;
    }

    @VisibleForTesting
    static boolean isClonal(final double copyNumber) {
        return Doubles.lessOrEqual(Doubles.absDistanceFromInteger(copyNumber), CLONAL_DISTANCE);
    }

}
