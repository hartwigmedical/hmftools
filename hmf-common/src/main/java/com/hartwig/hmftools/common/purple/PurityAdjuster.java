package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public abstract class PurityAdjuster {

    public static double impliedSamplePloidy(final double purity, final double normFactor) {
        return new PurityAdjusterTypicalChromosome(Gender.FEMALE, purity, normFactor).purityAdjustedCopyNumber("1", 1);
    }

    private final double purity;
    private final double normFactor;

    public PurityAdjuster(@NotNull final FittedPurity fittedPurity) {
        this(fittedPurity.purity(), fittedPurity.normFactor());
    }

    public PurityAdjuster(final double purity, final double normFactor) {
        this.purity = purity;
        this.normFactor = normFactor;
    }

    public double purity() {
        return purity;
    }

    public double normFactor() {
        return normFactor;
    }

    public abstract int germlineCopyNumber(@NotNull String contig);

    public double typicalRatio(@NotNull String contig) {
        return germlineCopyNumber(contig) / 2.0;
    }

    public double purityAdjustedCopyNumber(final String chromosomeName, final double ratio) {
        final double typicalRatio = typicalRatio(chromosomeName);
        return purityAdjustedCopyNumber(ratio, typicalRatio);
    }

    public double purityAdjustedCopyNumber(final double tumorRatio, final double normalRatio) {
        return Doubles.isZero(tumorRatio) ? 0 : 2 * normalRatio + 2 * (tumorRatio - normalRatio * normFactor) / purity / normFactor;
    }

    public double purityAdjustedVAF(@NotNull final String chromosome, final double copyNumber, final double observedFrequency) {
        int typicalCopyNumber = germlineCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, 0, copyNumber, observedFrequency);
    }

    public double purityAdjustedBAFSimple(final String chromosome, final double copyNumber, final double observedFrequency) {
        int typicalCopyNumber = germlineCopyNumber(chromosome);
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
        int typicalCopyNumber = germlineCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, 1, copyNumber, observedFrequency);
    }

    public double purityAdjustedVAFWithHomozygousNormal(@NotNull final String chromosome, final double copyNumber,
            final double observedFrequency) {
        int typicalCopyNumber = germlineCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, 2, copyNumber, observedFrequency);
    }
}
