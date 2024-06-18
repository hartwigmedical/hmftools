package com.hartwig.hmftools.wisp.purity.loh;

import static java.lang.Math.abs;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOD_MAX_ITERATIONS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOD_MIN_PROB_DIFF_PERC;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class LodCalcs
{
    public static double calcLimitOfDetection(double requiredProb, double meanCopyNumber, int lohFragments)
    {
        if(lohFragments == 0 || meanCopyNumber == 0)
            return 0;

        int iterations = 0;

        double expectedFragments = 0.5 * lohFragments;
        PoissonDistribution poisson = new PoissonDistribution(expectedFragments);

        double testValueLower = lohFragments < 10000 ? 0.005 : 10000.0 / lohFragments * 0.005;
        double testValueUpper = 1.0;

        double currentValue = (testValueLower + testValueUpper) * 0.5;

        int expectedValue = (int)round((1 - currentValue) / (currentValue * (meanCopyNumber - 2) + 2) * lohFragments);
        double currentProb = poisson.cumulativeProbability(expectedValue);
        double probDiff = 0;

        int lastUpperExpectedValue = 0;
        int lastLowerExpectedValue = 0;
        double lastUpperValue = 0;
        double lastLowerValue = 0;

        while(iterations < LOD_MAX_ITERATIONS)
        {
            probDiff = abs(requiredProb - currentProb);

            if(probDiff <= LOD_MIN_PROB_DIFF_PERC)
                break;

            // if prob is too low, need to lower the test value
            if(currentProb < requiredProb)
            {
                lastLowerExpectedValue = expectedValue;
                lastLowerValue = currentValue;
                testValueUpper = currentValue;
                currentValue = (currentValue + testValueLower) * 0.5;
            }
            else
            {
                lastUpperExpectedValue = expectedValue;
                lastUpperValue = currentValue;
                testValueLower = currentValue;
                currentValue =(currentValue + testValueUpper) * 0.5;
            }

            expectedValue = (int)round((1 - currentValue) / (currentValue * (meanCopyNumber - 2) + 2) * lohFragments);

            if(expectedValue == lastLowerExpectedValue || expectedValue == lastUpperExpectedValue)
            {
                currentValue = (lastLowerValue + lastUpperValue) * 0.5;
                break;
            }

            currentProb = poisson.cumulativeProbability(expectedValue);
            ++iterations;
        }

        if(iterations >= LOD_MAX_ITERATIONS)
        {
            CT_LOGGER.warn(format("max iterations reached: expectedFrags(%.0f) testValue(%.5f) prob(%.4f diff=%.4f)",
                    expectedFragments, currentValue, currentProb, probDiff));
        }

        return currentValue;
    }
}
