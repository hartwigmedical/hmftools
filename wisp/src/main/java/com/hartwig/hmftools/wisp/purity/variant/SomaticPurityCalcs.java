package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class SomaticPurityCalcs
{
    public static double estimatedPurityOld(final double tumorPurity, final double tumorVaf, final double sampleVaf, final double noiseRate)
    {
        return max(sampleVaf - noiseRate, 0) / tumorVaf * tumorPurity;
    }

    public static double estimatedPurity(final double sampleVaf, final double noiseRate, final FragmentTotals fragmentTotals)
    {
        double noiseAdjSampleVaf = max(sampleVaf - noiseRate, 0);
        double weightedVcn = fragmentTotals.weightedVariantCopyNumber();
        double weightedCn = fragmentTotals.weightedCopyNumber();
        return 2 * noiseAdjSampleVaf * (weightedVcn + 2 * noiseAdjSampleVaf - weightedCn * noiseAdjSampleVaf);
    }

    public static double estimatedProbability(final FragmentTotals fragmentTotals, double noiseRate)
    {
        int sampleFragments = fragmentTotals.sampleAdTotal();
        double expectedNoiseFragments = noiseRate * fragmentTotals.sampleDepthTotal();

        PoissonDistribution poisson = new PoissonDistribution(expectedNoiseFragments);

        if(sampleFragments == 0)
            return 1.0;

        double probability = 1 - poisson.cumulativeProbability(sampleFragments - 1);
        return probability;
    }

    public static double estimatedProbabilityNewNowOld(final FragmentTotals fragmentTotals, double noiseRate)
    {
        double weightedSampleDepth = fragmentTotals.weightedSampleDepth();
        double weightedSampleFrags = fragmentTotals.weightedSampleFrags();

        double effectiveDepth = weightedSampleDepth * fragmentTotals.variantCount();
        double effectiveAlleleFrags = weightedSampleFrags * fragmentTotals.variantCount();
        double expectedNoiseFragments = noiseRate * effectiveDepth;

        PoissonDistribution poisson = new PoissonDistribution(expectedNoiseFragments);

        if(effectiveAlleleFrags < 1)
        {
            return 1 - poisson.cumulativeProbability(0);
        }

        int lowerFrags = (int)floor(effectiveAlleleFrags - 1);
        int upperFrags = (int)ceil(effectiveAlleleFrags - 1);
        double probabilityLower = 1 - poisson.cumulativeProbability(lowerFrags);
        double probabilityUpper = 1 - poisson.cumulativeProbability(upperFrags);
        double upperWeight = effectiveAlleleFrags - floor(effectiveAlleleFrags);
        double lowerWeight = 1 - upperWeight;

        return lowerWeight * probabilityLower + upperWeight * probabilityUpper;
    }

    // OLD METHODS: likely deprecated
    public static double estimatedProbabilityOld(int alleleCount, double noise)
    {
        if(alleleCount <= noise || noise == 0)
            return 1;

        PoissonDistribution poisson = new PoissonDistribution(noise);
        double probability = 1 - poisson.cumulativeProbability(alleleCount - 1);
        return probability;
    }

    private static final int MAX_ITERATIONS = 20;
    private static final double MIN_PROB_DIFF_PERC = 0.001;

    public static double calcLimitOfDetection(final FragmentTotals fragmentTotals, final double noiseRate)
    {
        // first find the AF matches a specified limit of detection probability
        double sampleDepth = fragmentTotals.sampleDepthTotal();
        double expectedNoiseFrags = noiseRate * sampleDepth;

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

        while(iterations < MAX_ITERATIONS)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff <= MIN_PROB_DIFF_PERC)
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

            /*
            if(currentValue == testValueLower || currentValue == testValueUpper)
            {
                if(currentValue == testValueLower)
                    currentValue = newTestValue + 1;

                break;
            }
            */

            currentProb = 1 - poisson.cumulativeProbability(currentValue - 1);
            ++iterations;
        }

        if(iterations >= MAX_ITERATIONS)
        {
            CT_LOGGER.warn(format("max iterations reached: value(%.4f) test(%d) prob(%.4f diff=%.4f)",
                    expectedNoiseFrags, currentValue, currentProb, probDiff));
        }

        double af = currentValue / sampleDepth - noiseRate;

        double lod = 2 * af * (fragmentTotals.weightedVariantCopyNumber() + 2 * af - fragmentTotals.weightedCopyNumber() * af);

        return lod;
    }

    public static double calcLimitOfDetectionNewNowOld(final FragmentTotals fragmentTotals, final double noiseRate)
    {
        double effectiveFrags = fragmentTotals.weightedSampleDepth() * fragmentTotals.variantCount();
        double expectedNoiseFrags = noiseRate * effectiveFrags;

        int iterations = 0;
        double requiredProb = LOW_PROBABILITY;

        PoissonDistribution poisson = new PoissonDistribution(expectedNoiseFrags);

        int testValueLower = 1;
        int testValueUpper = (int)ceil(expectedNoiseFrags) * 10;
        int currentValue = testValueUpper;

        double currentProb = 1 - poisson.cumulativeProbability(currentValue - 1);
        double probDiff = 0;

        double lastProbUpper = 1.0;
        double lastProbLower = 0;

        while(iterations < MAX_ITERATIONS)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff <= MIN_PROB_DIFF_PERC)
            {
                if(currentProb < requiredProb)
                    ++currentValue;

                break;
            }

            // if prob is too low, need to lower the test value
            int newTestValue;
            if(currentProb < requiredProb)
            {
                lastProbLower = currentProb;
                testValueUpper = currentValue;
                newTestValue = (int)round((currentValue + testValueLower) * 0.5);
            }
            else
            {
                lastProbUpper = currentProb;
                testValueLower = currentValue;
                newTestValue = (int)round((currentValue + testValueUpper) * 0.5);
            }

            currentValue = newTestValue;

            if(currentValue == testValueLower || currentValue == testValueUpper)
            {
                if(currentValue == testValueLower)
                    currentValue = newTestValue + 1;

                break;
            }

            currentProb = 1 - poisson.cumulativeProbability(currentValue - 1);
            ++iterations;
        }

        if(iterations >= MAX_ITERATIONS)
        {
            CT_LOGGER.warn(format("max iterations reached: value(%.4f) test(%d) prob(%.4f diff=%.4f)",
                    expectedNoiseFrags, currentValue, currentProb, probDiff));
        }

        int lodFragmentsLower = currentValue;
        int lodFragmentsUpper = currentValue - 1;
        double lodTumorFractionLower = lodFragmentsLower / effectiveFrags - noiseRate;
        double lodTumorFractionUpper = lodFragmentsUpper / effectiveFrags - noiseRate;

        double weightedLod = ((requiredProb - lastProbUpper) * lodTumorFractionLower + (lastProbLower - requiredProb) * lodTumorFractionUpper)
                / (lastProbLower - lastProbUpper);

        return weightedLod;
    }
}
