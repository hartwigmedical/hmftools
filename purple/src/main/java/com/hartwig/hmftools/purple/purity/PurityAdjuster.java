package com.hartwig.hmftools.purple.purity;

import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class PurityAdjuster
{
    private final double mPurity;
    private final double mNormFactor;

    private final Map<String, Double> mObservedRatioMap;
    private final Gender mGender;

    public PurityAdjuster(final double purity, final double normFactor, final CobaltChromosomes cobaltChromosomes)
    {
        mPurity = purity;
        mNormFactor = normFactor;

        mObservedRatioMap = cobaltChromosomes.chromosomes().stream()
                .collect(Collectors.toMap(CobaltChromosome::contig, CobaltChromosome::actualRatio));

        mGender = cobaltChromosomes.gender();
    }

    public PurityAdjuster(final Gender gender, final double purity, final double normFactor)
    {
        mPurity = purity;
        mNormFactor = normFactor;
        mObservedRatioMap = null;
        mGender = gender;
    }

    public double purity()
    {
        return mPurity;
    }

    public double normFactor()
    {
        return mNormFactor;
    }

    public double germlineCopyNumber(final String chromosome)
    {
        return germlineRatio(chromosome) * 2;
    }

    public double purityAdjustedCopyNumber(final String chromosomeName, final double ratio)
    {
        final double typicalRatio = germlineRatio(chromosomeName);
        return purityAdjustedCopyNumber(ratio, typicalRatio);
    }

    public double purityAdjustedCopyNumber(final double tumorRatio, final double normalRatio)
    {
        return Doubles.isZero(tumorRatio) ? 0 : 2 * normalRatio + 2 * (tumorRatio - normalRatio * mNormFactor) / mPurity / mNormFactor;
    }

    public double purityAdjustedVAF(final String chromosome, final double copyNumber, final double observedFrequency)
    {
        return purityAdjustedVAF(chromosome, copyNumber, observedFrequency, false);
    }

    public double purityAdjustedVAF(
            final String chromosome, final double copyNumber, final double observedFrequency, boolean isGermlineHetDeletion)
    {
        double typicalCopyNumber = !isGermlineHetDeletion ? germlineCopyNumber(chromosome) : 1.0;
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

    private double germlineRatio(final String chromosome)
    {
        if(mObservedRatioMap != null)
            return mObservedRatioMap.getOrDefault(chromosome, 0d);

        // for testing only
        return HumanChromosome.fromString(chromosome).isDiploid(mGender) ? 1 : 0.5;
    }

    double purityAdjustedFrequency(final double normalCopyNumber, final double normalPloidy, final double tumorCopyNumber,
            final double observedFrequency)
    {
        return purityAdjustedPloidy(normalCopyNumber, normalPloidy, tumorCopyNumber, observedFrequency) / tumorCopyNumber;
    }

    public double purityAdjustedPloidy(final double normalCopyNumber, final double normalPloidy, final double tumorCopyNumber,
            final double observedFrequency)
    {
        double totalObservations = mPurity * tumorCopyNumber + normalCopyNumber * (1 - mPurity);
        double normalObservations = normalPloidy * (1 - mPurity);
        return (observedFrequency * totalObservations - normalObservations) / mPurity;
    }

    public double expectedFrequency(final double normalCopyNumber, final int normalPloidy, final double tumorCopyNumber,
            final double tumorPloidy)
    {
        if(Doubles.lessOrEqual(tumorCopyNumber, 0))
        {
            return 0;
        }

        double totalObservations = mPurity * tumorCopyNumber + normalCopyNumber * (1 - mPurity);
        double normalObservations = normalPloidy * (1 - mPurity);
        double tumorObservations = tumorPloidy * mPurity;

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

    public static double impliedPloidy(double averageRatio, double purity, double normFactor)
    {
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

}
