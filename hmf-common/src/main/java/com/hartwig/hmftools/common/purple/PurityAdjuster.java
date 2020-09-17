package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class PurityAdjuster {

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

    public int typicalCopyNumber(@NotNull String chromosome) {
        return HumanChromosome.fromString(chromosome).isDiploid(gender) ? 2 : 1;
    }

    public double typicalRatio(@NotNull String contig) {
        final Chromosome chromosome = HumanChromosome.fromString(contig);
        return  chromosome.isDiploid(gender) ? 1 : 0.5;
    }

    @SuppressWarnings("unused")
    public double impliedPloidy() {
        // Don't delete per request of Mr Jon Baber!!!
        return (1 - normFactor) / purity / normFactor * 2 + 2;
    }

    public double purityAdjustedCopyNumber(final String chromosomeName, final double ratio) {
        final double typicalRatio = typicalRatio(chromosomeName);
        return purityAdjustedCopyNumber(ratio, typicalRatio);
    }

    public double purityAdjustedCopyNumber(final double tumorRatio, final double normalRatio) {
        return Doubles.isZero(tumorRatio) ? 0 : 2 * normalRatio + 2 * (tumorRatio - normalRatio * normFactor) / purity / normFactor;
    }

    public double purityAdjustedVAF(@NotNull final String chromosome, final double copyNumber, final double observedFrequency) {
        int typicalCopyNumber = typicalCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, 0, copyNumber, observedFrequency);
    }

    public double purityAdjustedBAFSimple(final String chromosome, final double copyNumber, final double observedFrequency) {
        int typicalCopyNumber = typicalCopyNumber(chromosome);
        if (typicalCopyNumber < 2 || Doubles.lessOrEqual(copyNumber, 1)) {
            return 1;
        }
        return purityAdjustedFrequency(2, 1, copyNumber, observedFrequency);
    }

    double purityAdjustedFrequency(final int normalCopyNumber, final int normalPloidy, final double tumorCopyNumber,
            final double observedFrequency) {
        return purityAdjustedPloidy(normalCopyNumber, normalPloidy, tumorCopyNumber, observedFrequency) / tumorCopyNumber;
    }

    public double purityAdjustedPloidy(final int normalCopyNumber, final int normalPloidy, final double tumorCopyNumber,
            final double observedFrequency) {
        double totalObservations = purity * tumorCopyNumber + normalCopyNumber * (1 - purity);
        double normalObservations = normalPloidy * (1 - purity);
        return (observedFrequency * totalObservations - normalObservations) / purity;
    }

    public double expectedFrequency(final int normalCopyNumber, final int normalPloidy, final double tumorCopyNumber,
            final double tumorPloidy) {
        if (Doubles.lessOrEqual(tumorCopyNumber, 0)) {
            return 0;
        }

        double totalObservations = purity * tumorCopyNumber + normalCopyNumber * (1 - purity);
        double normalObservations = normalPloidy * (1 - purity);
        double tumorObservations = tumorPloidy * purity;

        return (normalObservations + tumorObservations) / totalObservations;
    }

    public double purityAdjustedVAFWithHeterozygousNormal(@NotNull final String chromosome, final double copyNumber,
            final double observedFrequency) {
        int typicalCopyNumber = typicalCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, 1, copyNumber, observedFrequency);
    }

    public double purityAdjustedVAFWithHomozygousNormal(@NotNull final String chromosome, final double copyNumber,
            final double observedFrequency) {
        int typicalCopyNumber = typicalCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, 2, copyNumber, observedFrequency);
    }
}
