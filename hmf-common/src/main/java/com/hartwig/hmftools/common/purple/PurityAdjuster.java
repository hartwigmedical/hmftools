package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;

import org.jetbrains.annotations.NotNull;

public class PurityAdjuster {

    private static final double AMBIGUOUS_BAF = 0.542;
    private static final double CLONAL_DISTANCE = 0.25;
    private static final double MAX_DEVIATION_ADJUSTMENT = 0.20;

    public static double impliedSamplePloidy(final double purity, final double normFactor) {
        return new PurityAdjuster(Gender.FEMALE, purity, normFactor).purityAdjustedCopyNumber("1", 1);
    }

    @NotNull
    private final Gender gender;
    private final double purity;
    private final double normFactor;

    public PurityAdjuster(@NotNull final Gender gender, @NotNull final FittedPurity fittedPurity) {
        this(gender, fittedPurity.purity(), fittedPurity.normFactor());
    }

    public PurityAdjuster(@NotNull final Gender gender, final double purity, final double normFactor) {
        this.gender = gender;
        this.purity = purity;
        this.normFactor = normFactor;
    }

    public double purity() {
        return purity;
    }

    public double normFactor() {
        return normFactor;
    }

    public double impliedPloidy() {
        return (1 - normFactor()) / purity() / normFactor() * 2 + 2;
    }

    public double purityAdjustedCopyNumber(final String chromosomeName, final double ratio) {
        final Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
        final double typicalRatio = chromosome.isHomologous(gender) ? 1 : 0.5;
        return purityAdjustedCopyNumber(ratio, typicalRatio);
    }

    public double purityAdjustedCopyNumber(final double tumorRatio, final double normalRatio) {
        return Doubles.isZero(tumorRatio) ? 0 : 2 * normalRatio + 2 * (tumorRatio - normalRatio * normFactor) / purity / normFactor;
    }

    @Deprecated
    public double purityAdjustedVAF(final double copyNumber, final double observedFrequency) {
        return purityAdjustedFrequency(copyNumber, observedFrequency, 0);
    }

    public double purityAdjustedVAF(final String chromosome, final double copyNumber, final double observedFrequency) {
        int normalCopyNumber = HumanChromosome.fromString(chromosome).isHomologous(gender) ? 2 : 1;
        return purityAdjustedFrequency(copyNumber, observedFrequency, normalCopyNumber, 0);
    }

    public double purityAdjustedBAF(final String chromosomeName, final double copyNumber, final double observedFrequency) {
        final Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
        if (!chromosome.isHomologous(gender) || Doubles.lessOrEqual(copyNumber, 0.5) || Doubles.isZero(observedFrequency)) {
            return 0;
        }

        if (Doubles.lessOrEqual(copyNumber, 1)) {
            return 1;
        }

        double typicalFrequency = 0.5;
        double rawAdjustedBaf = purityAdjustedFrequency(copyNumber, observedFrequency, typicalFrequency);

        int ploidy = (int) Math.round(copyNumber);
        if (Doubles.lessOrEqual(observedFrequency, AMBIGUOUS_BAF) && isClonal(copyNumber) && ploidy > 0) {
            int minBetaAllele = BAFUtils.minAlleleCount(ploidy);
            double modelBAF = BAFUtils.modelBAF(purity, ploidy, minBetaAllele);
            if (Doubles.lessThan(modelBAF, AMBIGUOUS_BAF)) {
                return (double) minBetaAllele / ploidy;
            }
        }

        return rawAdjustedBaf;
    }

    public double purityAdjustedMaxCopyNumberDeviation(double maxCopyNumberDeviation) {
        return maxCopyNumberDeviation * Math.max(1, MAX_DEVIATION_ADJUSTMENT / purity);
    }

    @VisibleForTesting
    static boolean isClonal(final double copyNumber) {
        return Doubles.lessOrEqual(Doubles.distanceFromInteger(copyNumber), CLONAL_DISTANCE);
    }

    @Deprecated
    private double purityAdjustedFrequency(final double copyNumber, final double observedFrequency, final double typicalFrequency) {
        return purityAdjustedFrequency(copyNumber, observedFrequency, 2, typicalFrequency);
    }

    private double purityAdjustedFrequency(final double tumorCopyNumber, final double observedFrequency, final int normalCopyNumber,
            final double normalFrequency) {
        assert (greaterThan(tumorCopyNumber, 0));
        assert (greaterThan(purity, 0));

        double normalTotalAllele = normalCopyNumber * (1 - purity);
        double normalBetaAllele = normalTotalAllele * normalFrequency;
        double tumorTotalAllele = tumorCopyNumber * purity;

        return (observedFrequency * (normalTotalAllele + tumorTotalAllele) - normalBetaAllele) / tumorCopyNumber / purity;
    }

}
