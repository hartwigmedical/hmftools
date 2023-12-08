package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_CHECK_THRESHOLD;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_REF_MIN_BIAS;
import static com.hartwig.hmftools.sage.SageConstants.STRAND_BIAS_REF_MIN_DEPTH;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class StrandBiasCalcs
{
    private final double[] mStrandBiasValues;
    private final double mMaxStrandBiasValue;

    private static final int MAX_AD_VS_PROB = 1000;
    private static final int MAX_ITERATIONS = 30;
    private static final double STRAND_BIAS_PROB = 0.0005;
    private static final double STRAND_BIAS_DIFF = 0.001;
    private static final double PROB_DIFF = 0.0001;
    private static final double EXPECTED_RATE = 0.5;
    private static final double INVALID_STRAND_BIAS = -1;

    public StrandBiasCalcs()
    {
        mStrandBiasValues = new double[MAX_AD_VS_PROB + 1];
        buildProbabilityCache();
        mMaxStrandBiasValue = mStrandBiasValues[mStrandBiasValues.length - 1];
    }

    public boolean isDepthBelowProbability(final StrandBiasData strandBiasDataAlt, final StrandBiasData strandBiasDataRef, boolean checkRef)
    {
        double strandBias = strandBiasDataAlt.bias();
        int depth = strandBiasDataAlt.depth();

        double minStrandBias = min(strandBias, 1 - strandBias);
        if(minStrandBias > STRAND_BIAS_CHECK_THRESHOLD)
            return false;

        if(checkRef)
        {
            // to use the ref there must be min depth observed
            if(strandBiasDataRef.forward() < STRAND_BIAS_REF_MIN_DEPTH || strandBiasDataRef.reverse() < STRAND_BIAS_REF_MIN_DEPTH)
                return false;

            double refBias = strandBiasDataRef.bias();

            double minRefBias = min(refBias, 1 - refBias);

            if(minRefBias < STRAND_BIAS_REF_MIN_BIAS && (refBias < 0.5) == (strandBias < 0.5))
                return false;
        }

        double requiredStrandBias = depth < mStrandBiasValues.length ? mStrandBiasValues[depth] : mMaxStrandBiasValue;

        if(requiredStrandBias == INVALID_STRAND_BIAS)
            return false;

        return Doubles.lessOrEqual(minStrandBias, requiredStrandBias);
    }

    private static final int INVALID_OBSERVED = -1;
    private static final double EPSILON = 1e-6;

    private void buildProbabilityCache()
    {
        int currentObserved = 0;

        for(int depth = 1; depth <= MAX_AD_VS_PROB; ++depth)
        {
            if(depth < 100)
            {
                int observed = findMinStrandBiasVsProb(depth, currentObserved);

                if(observed == INVALID_OBSERVED)
                {
                    mStrandBiasValues[depth] = INVALID_STRAND_BIAS;
                }
                else
                {
                    // the observed level is the first which exceeds the probability threshold, so set this to epsilon below
                    mStrandBiasValues[depth] = (observed - EPSILON) / depth;
                    currentObserved = observed;
                }
            }
            else
            {
                mStrandBiasValues[depth] = calcDepthVsStrandBias(depth);
            }
        }
    }

    private static int findMinStrandBiasVsProb(int depth, int startObserved)
    {
        BinomialDistribution distribution = new BinomialDistribution(depth, EXPECTED_RATE);

        int maxObserved = depth / 2;
        for(int observed = startObserved; observed <= maxObserved; ++observed)
        {
            double prob = distribution.cumulativeProbability(observed);

            if(prob > STRAND_BIAS_PROB)
            {
                if(observed == 0)
                    return INVALID_OBSERVED;
                else
                    return observed;
            }
        }

        return maxObserved;
    }

    private static double calcDepthVsStrandBias(int depth)
    {
        double lowerSb = 0.000;
        double upperSb = 0.5;
        double currentSb = (lowerSb + upperSb) * 0.5;

        BinomialDistribution distribution = new BinomialDistribution(depth, EXPECTED_RATE);

        int iterations = 0;
        int lastObserved = 0;
        while(iterations < MAX_ITERATIONS)
        {
            int observed = (int)round(depth * currentSb);

            if(observed == 0)
                return INVALID_STRAND_BIAS;

            if(lastObserved > 0 && observed == lastObserved)
                return currentSb;

            double prob = distribution.cumulativeProbability(observed);

            if(abs(prob - STRAND_BIAS_PROB) < PROB_DIFF)
                return currentSb;

            if(prob < STRAND_BIAS_PROB)
                lowerSb = currentSb;
            else
                upperSb = currentSb;

            if(abs(upperSb - lowerSb) < STRAND_BIAS_DIFF)
                return currentSb;

            currentSb = (lowerSb + upperSb) * 0.5;
            lastObserved = observed;
            ++iterations;
        }

        return currentSb;
    }
}
