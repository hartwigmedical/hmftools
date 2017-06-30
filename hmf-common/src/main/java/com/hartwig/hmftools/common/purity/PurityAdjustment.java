package com.hartwig.hmftools.common.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public class PurityAdjustment {

    public static double impliedSamplePloidy(final double purity, final double normFactor) {
        return new PurityAdjustment(Gender.FEMALE, purity, normFactor).purityAdjustedCopyNumber("1", 1);
    }

    @NotNull
    private final Gender gender;
    private final double purity;
    private final double normFactor;

    public PurityAdjustment(@NotNull final Gender gender, final double purity, final double normFactor) {
        this.gender = gender;
        this.purity = purity;
        this.normFactor = normFactor;
    }

    public double purityAdjustedCopyNumber(final String chromosome, final double ratio) {
        final double typicalCopyNumber = isMaleSexChromosome(chromosome) ? 1 : 2;
        return purityAdjustedCopyNumber(ratio, typicalCopyNumber);
    }

    private double purityAdjustedCopyNumber(final double ratio, final double typicalCopyNumber) {
        return Doubles.isZero(ratio) ? 0 : typicalCopyNumber + (2 * ratio - typicalCopyNumber * normFactor) / purity / normFactor;
    }

    public static double purityAdjustedVAF(final double purity, final double copyNumber, final double observedFrequency) {
        return purityAdjustedFrequency(purity, copyNumber, observedFrequency, 0);
    }

    public static double purityAdjustedBAF(final double purity, final double copyNumber, final double observedFrequency) {
        return purityAdjustedFrequency(purity, copyNumber, observedFrequency, 0.5);
    }

    private static double purityAdjustedFrequency(final double purity, final double copyNumber, final double observedFrequency,
            final double normalFrequency) {
        assert (greaterThan(copyNumber, 0));
        assert (greaterThan(purity, 0));

        double normalPloidy = 2 * (1 - purity);
        double tumorPloidy = copyNumber * purity;
        double normalAmount = 2 * (1 - purity) * normalFrequency;

        return (observedFrequency * (normalPloidy + tumorPloidy) - normalAmount) / copyNumber / purity;
    }

    private boolean isMaleSexChromosome(String chromosome) {
        return gender.equals(Gender.MALE) && (chromosome.equals("X") || chromosome.equals("Y"));
    }

}
