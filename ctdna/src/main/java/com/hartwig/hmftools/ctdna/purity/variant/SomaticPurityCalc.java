package com.hartwig.hmftools.ctdna.purity.variant;

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.hartwig.hmftools.common.stats.PoissonCalcs;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class SomaticPurityCalc
{
    public static final double LOW_PROBABILITY = 0.05;
    public static final double HIGH_PROBABILITY = 1 - LOW_PROBABILITY;
    private static final double OBSERVED_ZERO_LOW_MEAN = 3.0;

    public static FragmentCalcResult calc(double tumorPloidy, double tumorVaf, int totalCount, int alleleCount, double noise)
    {
        // ctDNA_TF = 2 * cfDNA_VAF / [ PLOIDY * ADJ_PRIMARY_VAF + cfDNA_VAF * ( 2 - PLOIDY)]
        // ADJ_PRIMARY_VAF= PRIMARY_VAF * [ PURITY*PLOIDY - 2*(1-PURITY)]/PURITY/PLOIDY

        double sampleVaf = totalCount > 0 ? max(alleleCount - noise, 0) / (double)totalCount : 0;
        double estimatedPurity = estimatedPurity(sampleVaf, tumorPloidy, tumorVaf);

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
            estimatedPurityLow = estimatedPurity(sampleVafLow, tumorPloidy, tumorVaf);

            double sampleVafHigh = max(highProbAlleleCount - noise, 0) / (double) totalCount;
            estimatedPurityHigh = estimatedPurity(sampleVafHigh, tumorPloidy, tumorVaf);
        }

        return new FragmentCalcResult(sampleVaf, estimatedPurity, probability, estimatedPurityLow, estimatedPurityHigh);
    }

    public static double calcPoissonNoiseValue(int alleleCount, double requiredProb)
    {
        if(alleleCount == 0)
            return requiredProb == LOW_PROBABILITY ? OBSERVED_ZERO_LOW_MEAN : 0;

        return PoissonCalcs.calcPoissonNoiseValue(alleleCount, requiredProb);
    }

    protected static double estimatedPurity(double sampleVaf, double tumorPloidy, double tumorVaf)
    {
        return max(min(2 * sampleVaf / (tumorPloidy * tumorVaf + sampleVaf * (2 - tumorPloidy)), 1), 0);
    }
}
