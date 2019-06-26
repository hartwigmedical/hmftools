package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public class ClonalityFactory {

    @NotNull
    private final Gender gender;
    private final double purity;
    private final double ploidyCutoff;

    private static final int TRIALS = 100_000;

    @NotNull
    public static ClonalityFactory fromPurityAdjuster(@NotNull PurityAdjuster purityAdjuster, double ploidyCutoff) {
        return new ClonalityFactory(purityAdjuster.gender(), purityAdjuster.purity(), ploidyCutoff);
    }

    public ClonalityFactory(@NotNull final Gender gender, final double purity, final double ploidyCutoff) {
        this.gender = gender;
        this.purity = purity;
        this.ploidyCutoff = ploidyCutoff;
    }

    @NotNull
    Clonality determineClonalityForVariant(@NotNull SomaticVariant variant) {
        double copyNumber = variant.adjustedCopyNumber();
        int typicalCopyNumber = HumanChromosome.fromString(variant.chromosome()).isDiploid(gender) ? 2 : 1;

        double monoploidProbability = purity / (purity * copyNumber + typicalCopyNumber * (1 - purity));
        double monoploidSamples = variant.totalReadCount() * monoploidProbability;
        double inconsistentSamples = Math.max(copyNumber, 0) * monoploidSamples;

        final BinomialDistribution inconsistentDistribution = new BinomialDistribution(TRIALS, Math.min(1, inconsistentSamples / TRIALS));
        if (variant.alleleReadCount() > inconsistentDistribution.inverseCumulativeProbability(0.999)) {
            return Clonality.INCONSISTENT;
        }

        final BinomialDistribution monoploidDistribution = new BinomialDistribution(TRIALS, Math.min(1, monoploidSamples / TRIALS));
        if (variant.alleleReadCount() < monoploidDistribution.inverseCumulativeProbability(0.001) || Doubles.lessThan(variant.ploidy(),
                ploidyCutoff)) {
            return Clonality.SUBCLONAL;
        }

        return Clonality.CLONAL;

    }
}
