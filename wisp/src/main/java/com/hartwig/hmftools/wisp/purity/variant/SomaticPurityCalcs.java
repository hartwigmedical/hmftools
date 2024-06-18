package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOD_MAX_ITERATIONS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOD_MIN_PROB_DIFF_PERC;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class SomaticPurityCalcs
{
    public static double estimatedPurity(final double sampleVaf, final double noiseRate, final FragmentTotals fragmentTotals)
    {
        double noiseAdjSampleVaf = max(sampleVaf - noiseRate, 0);
        double weightedVcn = fragmentTotals.weightedVariantCopyNumber();
        double weightedCn = fragmentTotals.weightedCopyNumber();
        return 2 * noiseAdjSampleVaf / (weightedVcn + 2 * noiseAdjSampleVaf - weightedCn * noiseAdjSampleVaf);
    }

    public static double estimatedProbability(final FragmentTotals fragmentTotals, double noiseRate)
    {
        int sampleFragments = fragmentTotals.sampleAdTotal();
        double expectedNoiseFragments = noiseRate * fragmentTotals.sampleDepthTotal();

        return estimatedProbability(sampleFragments, expectedNoiseFragments);
    }

    public static double estimatedProbability(int alleleFragments, double noiseFragments)
    {
        if(alleleFragments == 0 && noiseFragments == 0)
            return 1;

        PoissonDistribution poisson = new PoissonDistribution(noiseFragments);
        double probability = 1 - poisson.cumulativeProbability(alleleFragments - 1);
        return probability;
    }

    public static double calcLimitOfDetection(final FragmentTotals fragmentTotals, final double noiseRate)
    {
        // first find the AF matches a specified limit of detection probability
        double sampleDepth = fragmentTotals.sampleDepthTotal();
        double expectedNoiseFrags = noiseRate * sampleDepth;

        if(expectedNoiseFrags <= 0)
            return -1;

        int iterations = 0;
        double requiredProb = LOW_PROBABILITY;

        PoissonDistribution poisson = new PoissonDistribution(expectedNoiseFrags);

        int testValueLower = 1;
        int testValueUpper = (int)ceil(expectedNoiseFrags) * 10;
        int currentValue = testValueUpper;

        double currentProb = 1 - poisson.cumulativeProbability(currentValue - 1);
        double probDiff = 0;

        double lastProbUpperValue = 1.0;
        double lastProbLowerValue = 0;

        while(iterations < LOD_MAX_ITERATIONS)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff <= LOD_MIN_PROB_DIFF_PERC)
            {
                if(currentProb < requiredProb)
                    ++currentValue;

                break;
            }

            // if prob is too low, need to lower the test value
            int newTestValue;
            if(currentProb < requiredProb)
            {
                lastProbUpperValue = currentProb;
                testValueUpper = currentValue;
                newTestValue = (int)round((currentValue + testValueLower) * 0.5);
            }
            else
            {
                lastProbLowerValue = currentProb;
                testValueLower = currentValue;
                newTestValue = (int)round((currentValue + testValueUpper) * 0.5);
            }

            if(testValueLower >= testValueUpper - 1)
            {
                currentValue = abs(requiredProb - lastProbLowerValue) < abs(requiredProb - lastProbUpperValue) ? testValueLower : testValueUpper;
                break;
            }

            currentValue = newTestValue;

            currentProb = 1 - poisson.cumulativeProbability(currentValue - 1);
            ++iterations;
        }

        if(iterations >= LOD_MAX_ITERATIONS)
        {
            CT_LOGGER.warn(format("max iterations reached: value(%.4f) test(%d) prob(%.4f diff=%.4f)",
                    expectedNoiseFrags, currentValue, currentProb, probDiff));
        }

        double af = currentValue / sampleDepth - noiseRate;

        double lod = 2 * af / (fragmentTotals.weightedVariantCopyNumber() + 2 * af - fragmentTotals.weightedCopyNumber() * af);

        return lod;
    }
}
