package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
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

    public static double calcLimitOfDetection(final FragmentTotals fragmentTotals, final double tumorPurity, final double noiseRate)
    {
        double sampleDepthTotal = fragmentTotals.sampleAdjustedDepthTotal();
        double expectedNoiseFragments = noiseRate * sampleDepthTotal;
        double lodFragments = calcPoissonNoiseValue((int)round(expectedNoiseFragments), LOW_PROBABILITY);
        double lodSampleVaf = lodFragments / sampleDepthTotal;
        return estimatedPurity(tumorPurity, fragmentTotals.adjTumorVaf(), lodSampleVaf, noiseRate);
    }

    public static double estimatedProbability(int alleleCount, double noise)
    {
        if(alleleCount <= noise || noise == 0)
            return 1;

        PoissonDistribution poisson = new PoissonDistribution(noise);
        double probability = 1 - poisson.cumulativeProbability(alleleCount - 1);
        return probability;
    }
}
