package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.hartwig.hmftools.common.stats.PoissonCalcs;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class SomaticPurityCalc
{
    public static final double LOW_PROBABILITY = 0.05;
    public static final double HIGH_PROBABILITY = 1 - LOW_PROBABILITY;

    public static FragmentCalcResult calc(
            double tumorPloidy, double tumorVaf, int totalCount, int alleleCount, double noisePerMillion)
    {
        double noise = totalCount / 1000000.0 * noisePerMillion;

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
            double lowProbAlleleCount = PoissonCalcs.calcPoissonNoiseValue(alleleCount, HIGH_PROBABILITY);
            double highProbAlleleCount = PoissonCalcs.calcPoissonNoiseValue(alleleCount, LOW_PROBABILITY);

            double sampleVafLow = max(lowProbAlleleCount - noise, 0) / (double) totalCount;
            estimatedPurityLow = estimatedPurity(sampleVafLow, tumorPloidy, tumorVaf);

            double sampleVafHigh = max(highProbAlleleCount - noise, 0) / (double) totalCount;
            estimatedPurityHigh = estimatedPurity(sampleVafHigh, tumorPloidy, tumorVaf);
        }

        return new FragmentCalcResult(noise, sampleVaf, estimatedPurity, probability, estimatedPurityLow, estimatedPurityHigh);
    }

    private static double estimatedPurity(double sampleVaf, double tumorPloidy, double tumorVaf)
    {
        return max(min(2 * sampleVaf / (tumorPloidy * tumorVaf + sampleVaf * (2 - tumorPloidy)), 1), 0);
    }
}
