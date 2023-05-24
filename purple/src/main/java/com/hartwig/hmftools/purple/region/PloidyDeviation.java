package com.hartwig.hmftools.purple.region;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.Doubles.round;

import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import com.google.common.collect.Maps;
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

    private boolean mUseCache;
    private final Map<Double,Double> mCachedProbabilities;
    private AtomicLong mCalls;
    private long mCacheHits;

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

        mUseCache = false;
        mCachedProbabilities = Maps.newHashMap();
        mCacheHits = 0;
        mCalls = new AtomicLong(0);
    }

    public void setUseCache() { mUseCache = true; }

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

        if(!mUseCache)
            return calcDeviationProbability(deviationValue, ploidy);
        else
            return getDeviationProbability(deviationValue, ploidy);
    }

    private double calcDeviationProbability(double deviationValue, double ploidy)
    {
        mCalls.incrementAndGet();
        return 2 * mDistribution.cumulativeProbability(deviationValue) - 1 + max(-0.5 - ploidy, 0);
    }

    private synchronized double getDeviationProbability(double deviationValue, double ploidy)
    {
        double deviationValueRounded = round(deviationValue, 4);

        Double cachedProbability = mCachedProbabilities.get(deviationValueRounded);

        if(cachedProbability != null)
        {
            ++mCacheHits;
            return cachedProbability;
        }

        double deviationProbability = calcDeviationProbability(deviationValue, ploidy);
        mCachedProbabilities.put(deviationValueRounded, deviationProbability);
        return deviationProbability;
    }

    private double subMinAdditionalPenalty(double minPloidy, double ploidy)
    {
        return Math.min(mMajorAlleleSubMinAdditionalPenalty, max(0, -mMajorAlleleSubMinAdditionalPenalty * (ploidy - minPloidy)));
    }
}
