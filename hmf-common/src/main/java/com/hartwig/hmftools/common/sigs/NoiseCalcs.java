package com.hartwig.hmftools.common.sigs;

import static java.lang.Double.max;
import static java.lang.Math.abs;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import java.util.Map;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NoiseCalcs
{
    private static final Logger LOGGER = LogManager.getLogger(NoiseCalcs.class);

    public static final double POISSON_DEFAULT_PROBABILITY = 0.00001;
    public static final int POISSON_RANGE_DEFAULT_MIN_VALUE = 10;

    public static int calcRangeValue(final Map<Integer,Integer> rangeMap, int value)
    {
        return calcRangeValue(rangeMap, value);
    }

    public static int calcRangeValue(final Map<Integer,Integer> rangeMap, int value, double requiredProb, int minValue)
    {
        Integer rangeVal = rangeMap.get(value);
        if (rangeVal == null)
        {
            rangeVal = calcPoissonRangeGivenProb(value, requiredProb, minValue);
            rangeMap.put(value, rangeVal);
        }

        return rangeVal;
    }

    public static int calcPoissonRangeGivenProb(int value)
    {
        return calcPoissonRangeGivenProb(value, POISSON_DEFAULT_PROBABILITY, POISSON_RANGE_DEFAULT_MIN_VALUE);
    }

    public static int calcPoissonRangeGivenProb(int value, double requiredProb, int minValue)
    {
        // calculates a range around a given value where the probability of falling within that range is within the specified probability
        // the range is calculated using the lower value, which gives a more conservative estimation since the probabilities are not
        // symmetrical around the expected value - ie prob(X - value) < prob(x + value), but these differences diminish quickly for larger
        // expected values
        if(value <= minValue)
            return minValue;

        PoissonDistribution poisson = new PoissonDistribution(value);

        int maxIterations = 10;
        int iterations = 0;

        double initRange = 3.7 / sqrt(value); // works for requiredProb = 1e-4
        int testValue = (int)max(round(value * (1 - initRange)), 0);
        int testValueUpper = (int)max(round(value * (1 - initRange * 0.5)), 0);
        int testValueLower = (int)max(round(value * (1 - initRange * 2)), 0);

        double currentProb = poisson.cumulativeProbability(testValue);
        double probDiff = 0;

        while(iterations < maxIterations)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff < 0.1)
                break;

            // if prob is too high, need to lower the test value
            if(currentProb > requiredProb)
            {
                if(testValue <= testValueLower + 1)
                    break;

                testValueUpper = testValue;
                testValue = (int)round((testValue + testValueLower) * 0.5);
            }
            else
            {
                if(testValue >= testValueUpper - 1)
                    break;

                testValueLower = testValue;
                testValue = (int)round((testValue + testValueUpper) * 0.5);
            }

            currentProb = poisson.cumulativeProbability(testValue);
            ++iterations;
        }

        if(iterations >= maxIterations)
        {
            LOGGER.warn(String.format("max iterations reached: value(%d) test(%d) prob(%.4f diff=%.4f)",
                    value, testValue, currentProb, probDiff));
        }

        return value - testValue;
    }
}
