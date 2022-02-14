package com.hartwig.hmftools.sage.config;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;

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
    private static final double STRAND_BAIS_CHECK_THRESHOLD = 0.1;


    public StrandBiasCalcs()
    {
        mStrandBiasValues = new double[MAX_AD_VS_PROB + 1];
        buildProbabilityCache();
        mMaxStrandBiasValue = mStrandBiasValues[mStrandBiasValues.length - 1];
    }

    public boolean isDepthBelowProbability(double strandBias, int depth)
    {
        double minStrandBias = min(strandBias, 1 - strandBias);
        if(minStrandBias > STRAND_BAIS_CHECK_THRESHOLD)
            return false;

        double requiredStrandBias = depth < mStrandBiasValues.length ? mStrandBiasValues[depth] : mMaxStrandBiasValue;

        if(requiredStrandBias == INVALID_STRAND_BIAS)
            return false;

        return Doubles.lessOrEqual(minStrandBias, requiredStrandBias);
    }

    private void buildProbabilityCache()
    {
        int currentObserved = 0;

        for(int depth = 1; depth <= MAX_AD_VS_PROB; ++depth)
        {
            if(depth < 100)
            {
                int observed = findMinStrandBiasVsProb(depth, currentObserved);

                if(observed < 0)
                {
                    mStrandBiasValues[depth] = INVALID_STRAND_BIAS;
                }
                else
                {
                    mStrandBiasValues[depth] = observed / (double) depth;
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
                    return -1;
                else
                    return observed - 1;
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
