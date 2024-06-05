package com.hartwig.hmftools.common.stats;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class PoissonCalcs
{
    private static final Logger LOGGER = LogManager.getLogger(PoissonCalcs.class);

    private static final int DEFAULT_MAX_ITERATIONS = 20;
    private static final double MIN_PROB_DIFF_PERC_HIGH = 0.005;

    public static double calcPoissonNoiseValue(int observedCount, double requiredProb)
    {
        double defaultProbDiffThreshold = requiredProb < 0.01 ? requiredProb * 0.1 : MIN_PROB_DIFF_PERC_HIGH;

        return calcPoissonNoiseValue(observedCount, requiredProb, defaultProbDiffThreshold, DEFAULT_MAX_ITERATIONS);
    }

    public static double calcPoissonNoiseValue(int observedCount, double requiredProb, double probDiffThreshold, int maxIterations)
    {
        // calculate the mean of a Poisson distribution where the observed value would fall at the required probability
        if(observedCount < 0)
            return 0;

        if(observedCount == 0 && requiredProb > 0.5)
            return 0;

        // find the mean for a Poisson distribution where the observed fragment count is at the required probability level
        int iterations = 0;

        double currentValue = 0;
        double testValueUpper = 0;
        double testValueLower = 0;

        if(requiredProb > 0.5)
        {
            currentValue = observedCount * 0.5;
            testValueUpper = observedCount;
            testValueLower = currentValue * 0.01;
        }
        else
        {
            testValueUpper = findUpperStartValue(observedCount, requiredProb, maxIterations);
            testValueLower = observedCount;
            currentValue = (testValueUpper + testValueLower) * 0.5;
        }

        PoissonDistribution poisson = new PoissonDistribution(currentValue);

        double currentProb = poisson.cumulativeProbability(observedCount);
        double probDiff = 0;

        while(iterations < maxIterations)
        {
            probDiff = abs(requiredProb - currentProb);

            if(probDiff <= probDiffThreshold)
                return currentValue;

            // if prob is too low, need to lower the test value
            if(currentProb < requiredProb)
            {
                testValueUpper = currentValue;
                currentValue = (currentValue + testValueLower) * 0.5;
            }
            else
            {
                testValueLower = currentValue;
                currentValue = (currentValue + testValueUpper) * 0.5;
            }

            poisson = new PoissonDistribution(currentValue);
            currentProb = poisson.cumulativeProbability(observedCount);
            ++iterations;
        }

        if(iterations >= maxIterations)
        {
            LOGGER.warn(format("max iterations reached: value(%d) test(%.1f) prob(%.4f diff=%.4f)",
                    observedCount, currentValue, currentProb, probDiff));
        }

        return currentValue;
    }

    public static double findUpperStartValue(int observedCount, double requiredProb, int maxIterations)
    {
        double currentValue = max(1, observedCount * 2);

        for(int i = 0; i < maxIterations; ++i)
        {
            PoissonDistribution poisson = new PoissonDistribution(currentValue);

            double currentProb = poisson.cumulativeProbability(observedCount);

            if(currentProb < requiredProb)
                break;

            currentValue *= 2;
        }

        return currentValue;
    }
}
