package com.hartwig.hmftools.purple.segment;

import org.apache.commons.math3.distribution.BinomialDistribution;

public final class ExpectedBAF
{
    private static final double DEFAULT_PERCENT = 0.8;

    public static double expectedBAF(int averageDepth)
    {
        return expectedBAF(averageDepth, DEFAULT_PERCENT);
    }

    public static double expectedBAF(int averageDepth, double percent)
    {
        int minDepth = averageDepth / 2;
        int maxDepth = averageDepth * 3 / 2;
        double totalProbability = 0;

        for(int i = minDepth; i < maxDepth; i++)
        {
            final BinomialDistribution distribution = new BinomialDistribution(i, 0.5);
            double probability = 1d * distribution.inverseCumulativeProbability(percent) / i;
            totalProbability += probability;
        }

        return totalProbability / (maxDepth - minDepth);
    }
}
