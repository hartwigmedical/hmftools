package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_QUAL_NOISE_CUTOFF;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.MIN_QUAL_PER_AD;

import com.hartwig.hmftools.common.stats.PoissonCalcs;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class SomaticPurityCalcs
{
    public static final double LOW_PROBABILITY = 0.05;
    public static final double HIGH_PROBABILITY = 1 - LOW_PROBABILITY;
    private static final double OBSERVED_ZERO_LOW_MEAN = 3.0;

    public static double expectedNoise(final FragmentTotals fragmentTotals, final double noiseReadsPerMillion)
    {
        double qualPerAllele = fragmentTotals.qualPerAllele();

        double lowQualNoiseFactor = qualPerAllele < LOW_QUAL_NOISE_CUTOFF ?
                (LOW_QUAL_NOISE_CUTOFF - qualPerAllele) / (LOW_QUAL_NOISE_CUTOFF - MIN_QUAL_PER_AD) : 0;

        return fragmentTotals.sampleDepthTotal() / 1000000.0 * noiseReadsPerMillion + lowQualNoiseFactor;
    }

    public static double estimatedPurity(final double tumorPurity, final double tumorVaf, final double sampleVaf, final double sampleNoise)
    {
        return (sampleVaf - sampleNoise) / tumorVaf * tumorPurity;
    }

    public static double estimatedProbability(int alleleCount, double noise)
    {
        if(alleleCount <= noise || noise == 0)
            return 1;

        PoissonDistribution poisson = new PoissonDistribution(noise);
        double probability = 1 - poisson.cumulativeProbability(alleleCount - 1);
        return probability;
    }

    /*
    public static PurityCalcData calc(double tumorPloidy, double tumorVaf, int totalCount, int alleleCount, double noise)
    {
        // TFctDNA = [wVAFctDNA-ε]/ wVAFtissue / Puritytissue



        double sampleVaf = totalCount > 0 ? max(alleleCount - noise, 0) / (double)totalCount : 0;
        double estimatedPurity = estimatedPurityOld(sampleVaf, tumorPloidy, tumorVaf);

        double probability = 1;

        if(alleleCount > noise && noise > 0)
        {
            PoissonDistribution poisson = new PoissonDistribution(noise);
            probability = 1 - poisson.cumulativeProbability(alleleCount - 1);
        }

        double estimatedPurityLow = 0;
        double estimatedPurityHigh = 0;

        if(totalCount > 0)
        {
            double lowProbAlleleCount = calcPoissonNoiseValue(alleleCount, HIGH_PROBABILITY); // note reversal of prob thresholds is intended
            double highProbAlleleCount = calcPoissonNoiseValue(alleleCount, LOW_PROBABILITY);

            double sampleVafLow = max(lowProbAlleleCount - noise, 0) / (double) totalCount;
            estimatedPurityLow = estimatedPurityOld(sampleVafLow, tumorPloidy, tumorVaf);

            double sampleVafHigh = max(highProbAlleleCount - noise, 0) / (double) totalCount;
            estimatedPurityHigh = estimatedPurityOld(sampleVafHigh, tumorPloidy, tumorVaf);
        }

        return new PurityCalcData(sampleVaf, estimatedPurity, probability, estimatedPurityLow, estimatedPurityHigh);
    }

    public static double calcPoissonNoiseValue(int alleleCount, double requiredProb)
    {
        if(alleleCount == 0)
            return requiredProb == LOW_PROBABILITY ? OBSERVED_ZERO_LOW_MEAN : 0;

        return PoissonCalcs.calcPoissonNoiseValue(alleleCount, requiredProb);
    }

    protected static double estimatedPurityOld(double sampleVaf, double tumorPloidy, double tumorVaf)
    {
        return max(min(2 * sampleVaf / (tumorPloidy * tumorVaf + sampleVaf * (2 - tumorPloidy)), 1), 0);
    }
    */
}
