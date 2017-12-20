package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public class ClonalityFactory {
    private final PurityAdjuster purityAdjuster;
    private final double ploidyCutoff;

    private static final int TRIALS = 100_000;

    public ClonalityFactory(final PurityAdjuster purityAdjuster, final double ploidyCutoff) {
        this.purityAdjuster = purityAdjuster;
        this.ploidyCutoff = 0;
    }


    @NotNull
    public static Clonality fromPloidy(final double clonalCutoff, final double ploidy) {
        if (Doubles.isZero(clonalCutoff)) {
            return Clonality.UNKNOWN;
        } else if (Doubles.greaterOrEqual(ploidy, clonalCutoff)) {
            return Clonality.CLONAL;
        }
        return Clonality.SUBCLONAL;
    }

    @NotNull
    Clonality fromSample(@NotNull final PurityAdjustedSomaticVariant variant) {

        double purity = purityAdjuster.purity();
        double copyNumber = variant.adjustedCopyNumber();
        int typicalCopyNumber = purityAdjuster.typicalCopyNumber(variant.chromosome());

        double monoploidProbability = purity / (purity * copyNumber + typicalCopyNumber * (1 - purity));
        double monoploidSamples = variant.totalReadCount() * monoploidProbability;
        double inconsistentSamples = Math.max(copyNumber, 0) * monoploidSamples;

        final BinomialDistribution inconsistentDistribution = new BinomialDistribution(TRIALS, Math.min(1, inconsistentSamples / TRIALS));
        if (variant.alleleReadCount() > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
            return Clonality.INCONSISTENT;
        }

        final BinomialDistribution monoploidDistribution = new BinomialDistribution(TRIALS, Math.min(1, monoploidSamples / TRIALS));
        if (variant.alleleReadCount() < monoploidDistribution.inverseCumulativeProbability(0.01)
                && Doubles.greaterOrEqual(variant.ploidy(), ploidyCutoff)) {
            return Clonality.SUBCLONAL;
        }

        return Clonality.CLONAL;
    }

}
