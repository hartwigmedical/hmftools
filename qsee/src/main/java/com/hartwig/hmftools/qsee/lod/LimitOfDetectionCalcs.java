package com.hartwig.hmftools.qsee.lod;

import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.purple.Gender;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class LimitOfDetectionCalcs
{
    private static final int MIN_ALT_READS = 3;

    public static final double REF_SAMPLE_PURITY = 1.0;  // Assume no contamination from other samples

    public static double calcTheoreticalVafMutationOnOneAllele(double purity, HumanChromosome chromosome, Gender gender)
    {
        // Assume there is 1 mutation on 1 allele at a particular position.
        // With 0.5 purity and 2 alleles (diploid), the theoretical probability of observing the mutation on a read (i.e. theoretical VAF) is:
        //   0.5 purity * (1 mutation / 2 alleles) = 0.25

        int alleleCount = chromosome.isDiploid(gender) ? 2 : 1;
        return purity * (1.0 / alleleCount);
    }

    public static double calcMissedVariantLikelihood(Map<DepthRange,Integer> depthCounts, double theoreticalVaf)
    {
        // 'Likelihood' is a misnomer (should be 'probability'), but will leave as is because renaming is non-trivial

        int depthCountsTotal = depthCounts.values().stream().mapToInt(Integer::intValue).sum();
        if(depthCountsTotal == 0)
            return 0;

        double mvlhWeightedMean = 0;

        for(DepthRange depthRange : depthCounts.keySet())
        {
            int depthCount = depthCounts.get(depthRange);
            double weight = (double) depthCount / depthCountsTotal;

            if(weight == 0)
                continue;

            double mvlhAtDepth = 1;
            if(depthRange.midpoint() > 0 && theoreticalVaf > 0)
            {
                BinomialDistribution distribution = new BinomialDistribution(depthRange.midpoint(), theoreticalVaf);
                mvlhAtDepth = distribution.cumulativeProbability(MIN_ALT_READS - 1);
            }

            mvlhWeightedMean += mvlhAtDepth * weight;
        }

        return mvlhWeightedMean;
    }
}
