package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
import static com.hartwig.hmftools.common.stats.PoissonCalcs.findUpperStartValue;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_QUAL_NOISE_CUTOFF;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.MIN_QUAL_PER_AD;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class SomaticPurityCalcs
{
    public static double estimatedPurity(final double tumorPurity, final double tumorVaf, final double sampleVaf, final double noiseRate)
    {
        return max(sampleVaf - noiseRate, 0) / tumorVaf * tumorPurity;
    }

    public static double estimatedProbability(int alleleCount, double noise)
    {
        if(alleleCount <= noise || noise == 0)
            return 1;

        PoissonDistribution poisson = new PoissonDistribution(noise);
        double probability = 1 - poisson.cumulativeProbability(alleleCount - 1);
        return probability;
    }

    public static double calcLimitOfDetectionOld(final FragmentTotals fragmentTotals, final double tumorPurity, final double noiseRate)
    {
        double sampleDepthTotal = fragmentTotals.sampleAdjustedDepthTotal();
        double expectedNoiseFragments = noiseRate * sampleDepthTotal;
        double lodFragments = calcPoissonNoiseValue((int)round(expectedNoiseFragments), LOW_PROBABILITY);
        double lodSampleVaf = lodFragments / sampleDepthTotal;
        return estimatedPurity(tumorPurity, fragmentTotals.adjTumorVaf(), lodSampleVaf, noiseRate);
    }

    private static final int MAX_ITERATIONS = 20;
    private static final double MIN_PROB_DIFF_PERC = 0.001;

    public static double calcLimitOfDetection(final FragmentTotals fragmentTotals, final double noiseRate)
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
            CT_LOGGER.warn(format("max iterations reached: value(%.4f) test(%.1f) prob(%.4f diff=%.4f)",
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
