package com.hartwig.hmftools.purple.region;

import static java.lang.Math.max;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.NormalDistribution;

public class PloidyDeviation
{
    private final double mStandardDeviation;
    private final double mMinStandardDeviationPerPloidyPoint;
    private final NormalDistribution mDistribution;

    private final double mMajorAlleleSubOnePenaltyMultiplier;
    private final double mMajorAlleleSubMinAdditionalPenalty;
    private final double mMinDeviation;

    public PloidyDeviation(
            final double standardDeviation, final double minStandardDeviationPerPloidyPoint,
            final double majorAlleleSubOnePenaltyMultiplier, final double majorAlleleSubOneAdditionalPenalty, final double minDeviation)
    {
        mStandardDeviation = standardDeviation;
        mMinStandardDeviationPerPloidyPoint = minStandardDeviationPerPloidyPoint;
        mMajorAlleleSubOnePenaltyMultiplier = Math.abs(majorAlleleSubOnePenaltyMultiplier);
        mMajorAlleleSubMinAdditionalPenalty = Math.abs(majorAlleleSubOneAdditionalPenalty);
        mMinDeviation = minDeviation;
        mDistribution = new NormalDistribution();
    }

    public double majorAlleleDeviation(double purity, double normFactor, double ploidy)
    {
        final double majorAlleleDeviationMultiplier = Doubles.greaterThan(ploidy, 0) && Doubles.lessThan(ploidy, 1)
                ? max(1, mMajorAlleleSubOnePenaltyMultiplier * (1 - ploidy)) : 1;

        final double deviation = majorAlleleDeviationMultiplier * alleleDeviation(purity, normFactor, ploidy)
                + subMinAdditionalPenalty(1, ploidy);

        return max(deviation, mMinDeviation);
    }

    public double minorAlleleDeviation(double purity, double normFactor, double ploidy)
    {
        final double deviation = alleleDeviation(purity, normFactor, ploidy) + subMinAdditionalPenalty(0, ploidy);
        return max(deviation, mMinDeviation);
    }

    private double alleleDeviation(double purity, double normFactor, double ploidy)
    {
        final double ploidyDistanceFromInteger = Doubles.lessThan(ploidy, -0.5)
                ? 0.5 : Doubles.absDistanceFromInteger(ploidy);

        double standardDeviationsPerPloidy = max(mMinStandardDeviationPerPloidyPoint, purity * normFactor / 2 / mStandardDeviation);
        double deviationValue = ploidyDistanceFromInteger * standardDeviationsPerPloidy;

        return calcDeviationProbability(deviationValue, ploidy);
    }

    private double calcDeviationProbability(double deviationValue, double ploidy)
    {
        return 2 * mDistribution.cumulativeProbability(deviationValue) - 1 + max(-0.5 - ploidy, 0);
    }

    private double subMinAdditionalPenalty(double minPloidy, double ploidy)
    {
        return Math.min(mMajorAlleleSubMinAdditionalPenalty, max(0, -mMajorAlleleSubMinAdditionalPenalty * (ploidy - minPloidy)));
    }
}
