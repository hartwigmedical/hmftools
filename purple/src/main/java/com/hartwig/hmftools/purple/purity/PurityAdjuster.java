package com.hartwig.hmftools.purple.purity;

import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public abstract class PurityAdjuster
{
    public static double impliedPloidy(double averageRatio, double purity, double normFactor)
    {
        // Don't delete per request of Mr Jon Baber!!!
        return (averageRatio - normFactor) / purity / normFactor * 2 + 2;
    }

    public static double impliedNormFactor(double averageRatio, double purity, double tumorPloidy)
    {
        return 2 * averageRatio / (2 - 2 * purity + tumorPloidy * purity);
    }

    public static double impliedAverageRatio(double purity, double normFactor, double tumorPloidy)
    {
        return normFactor * (1 - purity + tumorPloidy / 2d * purity);
    }

    private final double purity;
    private final double normFactor;

    public PurityAdjuster(@NotNull final FittedPurity fittedPurity)
    {
        this(fittedPurity.purity(), fittedPurity.normFactor());
    }

    public PurityAdjuster(final double purity, final double normFactor)
    {
        this.purity = purity;
        this.normFactor = normFactor;
    }

    public double purity()
    {
        return purity;
    }

    public double normFactor()
    {
        return normFactor;
    }

    public double germlineCopyNumber(@NotNull String contig)
    {
        return germlineRatio(contig) * 2;
    }

    public abstract double germlineRatio(@NotNull String contig);

    public double purityAdjustedCopyNumber(final String chromosomeName, final double ratio)
    {
        final double typicalRatio = germlineRatio(chromosomeName);
        return purityAdjustedCopyNumber(ratio, typicalRatio);
    }

    public double purityAdjustedCopyNumber(final double tumorRatio, final double normalRatio)
    {
        return Doubles.isZero(tumorRatio) ? 0 : 2 * normalRatio + 2 * (tumorRatio - normalRatio * normFactor) / purity / normFactor;
    }

    public double purityAdjustedVAF(@NotNull final String chromosome, final double copyNumber, final double observedFrequency)
    {
        double typicalCopyNumber = germlineCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, 0, copyNumber, observedFrequency);
    }

    public double purityAdjustedBAFSimple(final String chromosome, final double copyNumber, final double observedFrequency)
    {
        double typicalCopyNumber = germlineCopyNumber(chromosome);
        if(typicalCopyNumber < 2 || Doubles.lessOrEqual(copyNumber, 1))
        {
            return 1;
        }
        return purityAdjustedFrequency(2, 1, copyNumber, observedFrequency);
    }

    double purityAdjustedFrequency(final double normalCopyNumber, final double normalPloidy, final double tumorCopyNumber,
            final double observedFrequency)
    {
        return purityAdjustedPloidy(normalCopyNumber, normalPloidy, tumorCopyNumber, observedFrequency) / tumorCopyNumber;
    }

    public double purityAdjustedPloidy(final double normalCopyNumber, final double normalPloidy, final double tumorCopyNumber,
            final double observedFrequency)
    {
        double totalObservations = purity * tumorCopyNumber + normalCopyNumber * (1 - purity);
        double normalObservations = normalPloidy * (1 - purity);
        return (observedFrequency * totalObservations - normalObservations) / purity;
    }

    public double expectedFrequency(final double normalCopyNumber, final int normalPloidy, final double tumorCopyNumber,
            final double tumorPloidy)
    {
        if(Doubles.lessOrEqual(tumorCopyNumber, 0))
        {
            return 0;
        }

        double totalObservations = purity * tumorCopyNumber + normalCopyNumber * (1 - purity);
        double normalObservations = normalPloidy * (1 - purity);
        double tumorObservations = tumorPloidy * purity;

        return (normalObservations + tumorObservations) / totalObservations;
    }

    public double purityAdjustedVAFWithHeterozygousNormal(@NotNull final String chromosome, final double copyNumber,
            final double observedFrequency)
    {
        double typicalCopyNumber = germlineCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, typicalCopyNumber / 2d, copyNumber, observedFrequency);
    }

    public double purityAdjustedVAFWithHomozygousNormal(@NotNull final String chromosome, final double copyNumber,
            final double observedFrequency)
    {
        double typicalCopyNumber = germlineCopyNumber(chromosome);
        return purityAdjustedFrequency(typicalCopyNumber, typicalCopyNumber, copyNumber, observedFrequency);
    }
}
